import pandas as pd
path = "GWAS_CP_all.txt"
df = pd.read_csv(path, sep="\t")

print("Number of Rows and Columns", df.shape)

print("Column Headers", list(df.columns))

print(df.head(4))