import pandas as pd
import numpy as np



df2 = pd.read_excel('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/core_regulatory_rules.xls', sheet_name='regulatory rules')

print(df2.columns)
df = df2[['Gene', 'Rule']]
df3 = df2[['bNum','Gene']]
#creating a dict from bNum to gene, and gene to bNum
gene_to_bnum = {}
bnum_to_gene = {}
for i in range(df3.shape[0]):
    gene_to_bnum[df3.iloc[i,1]] = df3.iloc[i,0]
    bnum_to_gene[df3.iloc[i,0]] = df3.iloc[i,1]
    


#saving the rules as a dictionary
rules = {}
for i in range(df.shape[0]):
    rules[df.iloc[i,0]] = df.iloc[i,1]


#saving rules as a json file
import json
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/core_regulatory_rules.json', 'w') as fp:
    json.dump(rules, fp)
#savinf gene to bnum as a json file
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/gene_to_bnum.json', 'w') as fp:
    json.dump(gene_to_bnum, fp)
#saving bnum to gene as a json file
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/bnum_to_gene.json', 'w') as fp:
    json.dump(bnum_to_gene, fp)

df = pd.read_excel('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/core_regulatory_rules.xls', sheet_name='metabolite_inputs')
df1 = pd.read_excel('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/core_regulatory_rules.xls', sheet_name='rxn_inputs')

met_list = df['metabolite inputs'] 
rxn_list = df1['reaction inputs']

#saving as csv files
met_list.to_csv('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/met_list1.csv')
rxn_list.to_csv('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/rxn_list.csv')