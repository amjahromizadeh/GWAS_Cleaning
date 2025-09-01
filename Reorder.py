import pandas as pd
path = "GWAS_CP_all.txt"
df = pd.read_csv(path, sep="\t", low_memory=False)

df = df.rename(columns={
    "MarkerName": "SNP",
    "EAF": "freq",
    "Beta": "b",
    "SE": "se",
    "Pval": "p"
})

df = df.drop(columns=["CHR", "POS"])

for col in ["freq", "b", "se", "p"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["n"] = 257841

df = df[["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]]

print("RnC =", df.shape)
print(df.head(4))

df.to_csv("author_formatted.tsv.gz", sep="\t", index=False, compression="gzip")