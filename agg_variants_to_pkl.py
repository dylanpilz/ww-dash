import pandas as pd

# Read in the data
df = pd.read_csv("app/frontend_data/aggregate_variants.csv")

df.to_pickle("app/frontend_data/aggregate_variants.pkl")