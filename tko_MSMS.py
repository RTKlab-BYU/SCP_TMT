
# Import packages
import pandas as pd
import plotly.express as px

important_proteins = ["MET6","HIS4","URA2"]
# file_name = "peptides 2.txt"
file_name = "msms 1.txt"
run_filter = "Quan"
colors = ["blue","blue","blue","red","red","red","green","green","green","gray","gray"]
# run_filter = "id_"


# tko_channel_groups = {"MET6 KO":[1,2,3],
#                       "HIS4 KO":[4,5,6],
#                       "URA2 KO":[7,8,9],
#                       "Blank":[10,11]}

# filter rows
data = pd.read_table(file_name,sep="\t")
data = data[data["Gene Names"].isin(important_proteins)]
data = data.loc[data["Raw file"].str.contains(run_filter)]

#filter columns
prot_abundance = data.filter(regex='Reporter intensity |Gene Names|Raw file')
prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity corrected ")]
prot_abundance = prot_abundance.loc[:,~prot_abundance.columns.str.contains("Reporter intensity count ")]
# print(data["Gene names"])
# print(prot_abundance)

#melt
numeric_cols =  prot_abundance.columns
numeric_cols = numeric_cols[numeric_cols != "Gene Names"]
prot_abundance = pd.melt(prot_abundance,id_vars=["Gene Names","Raw file"], value_vars = numeric_cols,var_name="Channel",value_name="intensity")
# prot_abundance = prot_abundance.loc[prot_abundance["intensity"]!=0]
print(prot_abundance)
prot_abundance["Channel"] = prot_abundance["Channel"].str[19:].astype(int)
# print(prot_abundance)

#summarize
# prot_abundance = prot_abundance.groupby(by =["Gene Names","Channel"]).agg({"intensity":[("intensity","sum"),("stdev","std")]})
prot_abundance = prot_abundance.groupby(by =["Gene Names","Channel","Raw file"]).agg({"intensity":[("intensity","sum")]})
# prot_abundance = prot_abundance.groupby(by =["Gene Names","Channel"]).agg({"intensity":[("intensity","sum")]})
prot_abundance.columns= prot_abundance.columns.get_level_values(1)
prot_abundance = prot_abundance.reset_index()
prot_abundance = prot_abundance.groupby(by =["Gene Names","Channel"]).agg({"intensity":[("intensity","mean"),("stdev","std")]})
prot_abundance.columns= prot_abundance.columns.get_level_values(1)
prot_abundance = prot_abundance.reset_index()
print(prot_abundance)


#plot
for each_protein in important_proteins:
    fig = px.bar(prot_abundance[prot_abundance["Gene Names"]==each_protein],
                #   x="Channel", y = "intensity",color=colors,title=each_protein)  
                  x="Channel", y = "intensity",error_y="stdev",color=colors,title=each_protein)  
    fig.show()