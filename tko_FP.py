import pandas as pd
import plotly.express as px
important_proteins = ["MET6","HIS4","URA2"]
file_name = "ratio_protein_MD.tsv"
num_channels = 11
# run_filter = "id_"


# tko_channel_groups = {"MET6 KO":[1,2,3],
#                       "HIS4 KO":[4,5,6],
#                       "URA2 KO":[7,8,9],
#                       "Blank":[10,11]}

data = pd.read_table(file_name,sep="\t")
data = data[data["Gene"].isin(important_proteins)]
prot_abundance = data.filter(regex='sample|Gene')
# print(data["Gene"])
print(prot_abundance)
numeric_cols =  prot_abundance.columns
numeric_cols = numeric_cols[numeric_cols != "Gene"]
# prot_abundance.columns = prot_abundance.columns.rename({prot_abundance.columns,prot_abundance.columns.str.replace("\\.[d]+","",regex=True)})
prot_abundance = pd.melt(prot_abundance,id_vars="Gene", value_vars = numeric_cols,var_name="reporter",value_name="intensity")
# prot_abundance["reporter"] = prot_abundance["reporter"].str[19:21].astype(int)
# print(prot_abundance)
prot_abundance = prot_abundance.groupby(by =["Gene","reporter"]).agg({"intensity":[("intensity","mean"),("stdev","std")]})
prot_abundance.columns= prot_abundance.columns.get_level_values(1)
prot_abundance = prot_abundance.reset_index()
print(prot_abundance)

for each_protein in important_proteins:
    fig = px.bar(prot_abundance[prot_abundance["Gene"]==each_protein],
                  x="reporter", y = "intensity", error_y="stdev",title=each_protein)  
    fig.show()