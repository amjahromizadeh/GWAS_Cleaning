import pandas as pd
from pathlib import Path

# 1) File paths (adjust if yours live elsewhere)
panel_path = "ref_panel/integrated_call_samples_v3.20130502.ALL.panel"
related_path = "ref_panel/20140625_related_individuals.txt"  # optional

# 2) Read as "any whitespace" separated (tabs/spaces)
#    - use a RAW string for the regex: r"\s+"
#    - engine="python" supports regex separators cleanly
panel = pd.read_csv(panel_path, sep=r"\s+", engine="python")

# 3) Normalise header names to lowercase
panel.columns = [c.lower() for c in panel.columns]

# 4) Find the super-population column (support both spellings)
if "super_pop" in panel.columns:
    sp_col = "super_pop"
elif "super_population" in panel.columns:
    sp_col = "super_population"
else:
    raise KeyError(f"Can't find a super-pop column. Columns are: {list(panel.columns)}")

# 5) Keep EUR super-population
eur = panel[panel[sp_col].astype(str).str.upper() == "EUR"].copy()

# 6) (Optional) remove related individuals if the file exists
try:
    # Read first column only (sample IDs), whitespace-separated
    related = pd.read_csv(related_path, sep=r"\s+", engine="python",
                          header=None, usecols=[0], names=["sample"])
    related["sample"] = related["sample"].astype(str).str.strip()
    eur = eur[~eur["sample"].isin(set(related["sample"]))]
    print(f"Excluded related individuals listed in: {related_path}")
except FileNotFoundError:
    print("No related list found; keeping all EUR samples.")

# 7) Write one sample ID per line (exactly as in VCF headers)
out_dir = Path("samples")
out_dir.mkdir(exist_ok=True)
eur_path = out_dir / "EUR.unrelated.samples"
eur["sample"].to_csv(eur_path, index=False, header=False)

# 8) Quick checks
print("Wrote:", eur_path)
print("# EUR samples:", eur.shape[0])
print("First few:\n", eur["sample"].head(5).to_string(index=False))
