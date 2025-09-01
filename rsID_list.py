import pandas as pd

path = "author_formatted.tsv.gz"

rsids = pd.read_csv(path, sep="\t", usecols=["SNP"])

rsids["SNP"] = rsids["SNP"].astype(str).str.strip().str.lower()

rsids = rsids[rsids["SNP"].str.match(r"^rs\d+$", na=False)]

rsids = rsids.drop_duplicates()

out_path = "rsids.txt"
rsids["SNP"].to_csv("rsids.txt", index=False, header=False)

mask = rsids["SNP"].str.match(r"^rs\d+$", na=False)
kept = rsids[mask]
dropped = rsids[~mask]


print("Total Unique rsIDs:", len(rsids))
print(rsids.head(5).to_string(index=False))
print("Dropped: ", dropped.head(5).to_string(index=False))