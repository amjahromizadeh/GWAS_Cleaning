import pandas as pd

# 1) Point to inputs
gwas_path = "author_formatted.tsv.gz"   # your formatted GWAS: SNP, A1, A2, freq, b, se, p, n
ref_path  = "ref_panel.tsv.gz"          # your reference table: CHROM, POS, ID, REF, ALT, AF_EUR

# 2) Read both tables
gwas = pd.read_csv(gwas_path, sep="\t")
ref  = pd.read_csv(ref_path,  sep="\t")

# 3) Normalize key columns so they match exactly during merge
gwas["SNP"] = gwas["SNP"].astype(str).str.strip().str.lower()
ref["ID"]   = ref["ID"].astype(str).str.strip().str.lower()

# 4) (Safety) drop any duplicate IDs in the reference to avoid row duplication on merge
#    - keep='first' is fine because the ref was constructed as biallelic SNPs
dup_n = ref.duplicated(subset=["ID"]).sum()
if dup_n:
    print(f"Warning: reference has {dup_n} duplicate rsIDs; keeping the first occurrence.")
    ref = ref.drop_duplicates(subset=["ID"], keep="first")

# 5) Keep only the columns we want from the reference
ref_small = ref[["ID", "CHROM", "POS", "REF", "ALT", "AF_EUR"]].copy()

# 6) Left-merge: keep ALL GWAS rows; add reference columns when rsID matches
merged = gwas.merge(ref_small, left_on="SNP", right_on="ID", how="left")

# 7) We no longer need the helper 'ID' column after the merge
merged = merged.drop(columns=["ID"])

# 8) Arrange columns in a tidy order
merged = merged[["SNP", "A1", "A2", "freq", "b", "se", "p", "n",
                 "CHROM", "POS", "REF", "ALT", "AF_EUR"]]

# 9) Save the merged table
out_path = "gwas_plus_ref.tsv.gz"
merged.to_csv(out_path, sep="\t", index=False, compression="gzip")

# 10) Quick sanity stats + write unmatched rsIDs for inspection
total = len(merged)
matched = merged["REF"].notna().sum()  # REF present implies a match to the reference
print(f"Matched to reference: {matched}/{total} ({matched/total:.1%})")
unmatched = merged.loc[merged["REF"].isna(), "SNP"].dropna().drop_duplicates()
unmatched.to_csv("unmatched_rsids.txt", index=False, header=False)
print(f"Unmatched rsIDs written to unmatched_rsids.txt: {len(unmatched)}")
print("Wrote:", out_path)
