import pandas as pd
from pathlib import Path

core = pd.read_csv("ref_panel/ref/ref_core.tsv", sep="\t")  # has CHROM, POS, ID, REF, ALT

# ensure types are clean
core["CHROM"] = core["CHROM"].astype(str)
core["POS"]   = core["POS"].astype(int)

# autosomes only (1..22)
auto = core[core["CHROM"].isin([str(i) for i in range(1, 23)])].copy()

outdir = Path("regions"); outdir.mkdir(exist_ok=True)

# For each chromosome, write a TAB-separated file with two columns: CHROM POS
for chrom, sub in auto.groupby("CHROM"):
    sub = sub.sort_values("POS")
    # build a 2-column DataFrame: CHROM, POS
    reg = pd.DataFrame({"CHROM": chrom, "POS": sub["POS"]})
    reg.to_csv(outdir / f"chr{chrom}.regions.txt",
               sep="\t", header=False, index=False)

print("Wrote per-chromosome region files (TAB-separated) to ./regions/")
