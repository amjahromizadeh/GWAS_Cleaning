import pandas as pd
path = "ref_panel.tsv.gz"

df = pd.read_csv(path, sep="\t", low_memory=False)

print("Columns: ", list(df.columns))
print("Shape: ", df.shape)
print(df.head(4))