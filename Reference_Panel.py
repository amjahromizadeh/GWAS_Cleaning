import time                         # lets us pause if the server asks us to slow down
import math                         # for ceiling when we compute number of batches
import json                         # to parse JSON if we need to print/debug
import requests                     # HTTP library – to call the REST API
import pandas as pd                 # table handling
from collections import defaultdict # handy for building dicts with default values

# 0) Inputs and small knobs you can tweak
gwas_path = "author_formatted.tsv.gz"  # your GWAS with the 'SNP' column
batch_size = 200                       # Ensembl POST limit (documented)
server = "https://grch37.rest.ensembl.org"
endpoint = "/variation/homo_sapiens"   # POST endpoint for many IDs
params = {"pops": 1}                   # ask for population allele frequencies
headers = {"Content-Type": "application/json", "Accept": "application/json"}

# 1) Read your GWAS and collect unique rsIDs
gwas = pd.read_csv(gwas_path, sep="\t", low_memory=False)
rsids = pd.unique(gwas["SNP"].astype(str))   # unique SNP IDs as strings
print("Unique rsIDs:", len(rsids))

# 2) Helper to parse one variant’s JSON into a flat record
def parse_variant(rsid, payload):
    """
    rsid: e.g., 'rs2352974'
    payload: the JSON object returned for that rsID
    returns: dict with SNP, chr, pos, ref, alt, strand, AF_1KG_ALL, AF_1KG_EUR
    """
    rec = {"SNP": rsid,
           "chr": None, "pos": None, "ref": None, "alt": None, "strand": None,
           "AF_1KG_ALL": None, "AF_1KG_EUR": None}

    # --- (A) pick the GRCh37 mapping and allele_string ---
    # The 'mappings' list contains positions per assembly; we want GRCh37.
    # Each mapping has: assembly_name, seq_region_name (chr), start, end, strand, allele_string.
    mappings = payload.get("mappings", [])
    grch37_maps = [m for m in mappings if m.get("assembly_name") == "GRCh37"]

    if grch37_maps:
        m = grch37_maps[0]  # take the first GRCh37 mapping
        rec["chr"] = str(m.get("seq_region_name"))
        rec["pos"] = int(m.get("start")) if m.get("start") is not None else None
        rec["strand"] = int(m.get("strand")) if m.get("strand") is not None else None

        # allele_string is typically "REF/ALT" (reference allele first, per Ensembl’s format docs)
        allele_string = m.get("allele_string")
        if allele_string and "/" in allele_string:
            left, right = allele_string.split("/", 1)
            rec["ref"] = left.upper()
            rec["alt"] = right.upper()

    # --- (B) scan population frequencies for 1000 Genomes Phase 3 ---
    # The 'populations' list holds entries with 'population', 'allele', 'frequency'.
    # We’ll pull overall ALL and EUR AFs where available.
    pops = payload.get("populations", [])
    af_all = None
    af_eur = None
    # We want frequency of the ALT allele by default (common for references).
    # But to make it simple & robust, we’ll store the frequency of whatever equals rec["alt"].
    for p in pops:
        pop_name = str(p.get("population", ""))
        allele = str(p.get("allele", "")).upper()
        freq = p.get("frequency", None)

        # 1000 Genomes Phase 3 labels look like "1000GENOMES:phase_3:ALL" or "1000GENOMES:phase_3:EUR"
        if pop_name.endswith(":ALL") and allele == rec["alt"]:
            af_all = freq if freq is not None else af_all
        if pop_name.endswith(":EUR") and allele == rec["alt"]:
            af_eur = freq if freq is not None else af_eur

    rec["AF_1KG_ALL"] = af_all
    rec["AF_1KG_EUR"] = af_eur
    return rec

# 3) Loop over rsIDs in batches and call the API
records = []                         # we’ll append one parsed dict per rsID
n_batches = math.ceil(len(rsids) / batch_size)

for i in range(n_batches):
    chunk = rsids[i*batch_size : (i+1)*batch_size]
    body = {"ids": [str(x) for x in chunk]}  # POST body: a list of rsIDs

    # Make the request; we add params (?pops=1) to include population allele freqs
    resp = requests.post(
        server + endpoint,
        params=params,
        headers=headers,
        data=json.dumps(body),
        timeout=60
    )

    # Basic politeness + backoff on HTTP 429 or if headers suggest we do so
    # (If rate-limited, Ensembl may set 'Retry-After'; if present, we sleep that many seconds.)
    if resp.status_code == 429:
        retry_after = int(resp.headers.get("Retry-After", "1"))
        time.sleep(retry_after)
        resp = requests.post(server + endpoint, params=params, headers=headers, data=json.dumps(body), timeout=60)

    resp.raise_for_status()
    payload = resp.json()  # dict mapping each id -> its JSON object (or error)

    # Parse each rsID in the chunk
    for rsid in chunk:
        item = payload.get(rsid)
        if isinstance(item, dict) and not item.get("error"):
            rec = parse_variant(rsid, item)
            records.append(rec)
        else:
            # Keep a placeholder so we know which ones failed (chr/ref/alt will be None)
            records.append({"SNP": rsid, "chr": None, "pos": None, "ref": None, "alt": None,
                            "strand": None, "AF_1KG_ALL": None, "AF_1KG_EUR": None})

    # Friendly progress print
    print(f"Batch {i+1}/{n_batches} done")

# 4) Turn into a DataFrame and save
refdf = pd.DataFrame.from_records(records)
print(refdf.head(5))
refdf.to_csv("ensembl_ref_grch37.tsv.gz", sep="\t", index=False, compression="gzip")
print("Saved: ensembl_ref_grch37.tsv.gz")
