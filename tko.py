import pandas as pd
import plotly.express as px
important_proteins = ["MET6","PFK2","URA2"]
file_name = "proteinGroups 2.txt"
run_filter = "Quan"

# tko_channel_groups = {"MET6 KO":[1,2,3],
#                       "PFK2 KO":[4,5,6],
#                       "URA2 KO":[7,8,9],
#                       "Blank":[10,11]}

data = pd.read_table(file_name,sep="\t")
data = data[data["Gene names"].isin(important_proteins)]
prot_abundance = data.filter(regex='Reporter intensity |Gene names')
prot_abundance = prot_abundance.filter(regex=run_filter+"|Gene names")
prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity corrected ")]
prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity count ")]
# print(data["Gene names"])
# print(prot_abundance)
numeric_cols =  prot_abundance.columns
numeric_cols = numeric_cols[numeric_cols != "Gene names"]
prot_abundance = pd.melt(prot_abundance,id_vars="Gene names", value_vars = numeric_cols,var_name="Run and reporter",value_name="intensity")
prot_abundance["Run and reporter"] = prot_abundance["Run and reporter"].str[19:21].astype(int)
# print(prot_abundance)
prot_abundance = prot_abundance.groupby(by =["Gene names","Run and reporter"]).agg({"intensity":[("intensity","mean"),("stdev","std")]})
prot_abundance.columns= prot_abundance.columns.get_level_values(1)
prot_abundance = prot_abundance.reset_index()
print(prot_abundance)

for each_protein in important_proteins:
    fig = px.bar(prot_abundance[prot_abundance["Gene names"]==each_protein],
                  x="Run and reporter", y = "intensity", error_y="stdev",title=each_protein)  
    fig.show()