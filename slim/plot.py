import plotly.express as px
import pandas as pd

df = pd.read_csv("tsv_output/dz1_mean_by_psi11_qtl.csv")
# df.drop(columns=["Unnamed: 0"], inplace=True)
df.set_index("qtl_count", inplace=True)
fig = px.imshow(df)
fig.show()