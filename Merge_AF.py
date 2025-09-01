import pandas as pd

# --- 1) point to inputs
core_path = "ref_panel/ref/ref_core.tsv"                # CHROM POS ID REF ALT
eur_path  = "ref_af/ref_eur_af.tsv"      # CHROM POS ID REF ALT AF_EUR

# --- 2) read files (only needed columns from AF file)
core = pd.read_csv(core_path, sep="\t")
eur  = pd.read_csv(eur_path,  sep="\t", usecols=["CHROM","POS","ID","REF","ALT","AF_EUR"])

# --- 3) normalize key columns so they match exactly
for df in (core, eur):
    df["CHROM"] = df["CHROM"].astype(str)          # e.g. "1".."22"
    df["POS"]   = df["POS"].astype(int)            # 1-based position
    df["ID"]    = df["ID"].astype(str).str.lower() # "rs123..."
    df["REF"]   = df["REF"].astype(str).str.upper()
    df["ALT"]   = df["ALT"].astype(str).str.upper()

# --- 4) (optional) sanity: drop dup rows in AF table on the join keys
eur = eur.drop_duplicates(subset=["CHROM","POS","ID","REF","ALT"])

# --- 5) left-merge: keep all core rows; add AF_EUR when there is a match
ref_panel = core.merge(
    eur,
    on=["CHROM","POS","ID","REF","ALT"],
    how="left"
)

# --- 6) save and quick peek
out_path = "ref_panel.tsv.gz"
ref_panel.to_csv(out_path, sep="\t", index=False, compression="gzip")

print("Rows:", len(ref_panel))
print("Columns:", list(ref_panel.columns))
print("AF_EUR non-missing:", ref_panel["AF_EUR"].notna().sum())
print(ref_panel.head(5).to_string(index=False))
print("Wrote:", out_path)
