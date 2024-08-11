import cobra 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr

# Load the model
model = cobra.io.read_sbml_model("/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/ecoli_core.xml")

#initalise the bounds here!!!!





# Load other data

#read from json as dictionary
import json
genes_to_rule = {}        ##################################################
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/core_regulatory_rules.json', 'r') as fp:
    genes_to_rule = json.load(fp)

for key, value in genes_to_rule.items():

    #converts values strings to exprs
    genes_to_rule[key] = parse_expr(value.replace("NOT", "~").replace("AND", "&").replace("OR", "|").replace("[", "_").replace("]", "").replace("-", "_"))

#read from csv as list
met_df = pd.read_csv('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/met_list.csv')
rxn_df = pd.read_csv('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/rxn_list.csv')

rxn_list = rxn_df['reaction inputs'].tolist() ##################################################

met_list = met_df['metabolite inputs'].tolist() ##################################################

symbol_translation_backwards = {}  ##################################################
symbol_translation_forwards = {}   ##################################################

for i in range(len(met_df['translated'])):
    symbol_translation_backwards[met_df["translated"][i]] = met_df["metabolite inputs"][i]
    symbol_translation_forwards[met_df["metabolite inputs"][i]] = met_df["translated"][i]

genes_to_var = {}  ##################################################
for x in genes_to_rule.keys():
    genes_to_var[x] = 1

#initiliazing genes to vars
genes_to_var['aceA'] = 0
genes_to_var['aceB'] = 1
genes_to_var['aceE'] = 0
#....

#making a dictionary of metabolites to inital concentrations
initial_concs_met = {}  ##################################################
for met in met_list:
    initial_concs_met[met] = 0

#open bnum to gene and gene to bnum
gene_to_bnum = {}  ##################################################
bnum_to_gene = {}  ##################################################
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/gene_to_bnum.json', 'r') as fp:
    gene_to_bnum = json.load(fp)
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/bnum_to_gene.json', 'r') as fp:
    bnum_to_gene = json.load(fp)



#making a dictionary of reactions in rxnlist to irules
rxns_to_rule = {}
i=0
for rxn in rxn_list:
    print(i)
    i+=1
    rx = model.reactions.get_by_id(rxn).gene_reaction_rule
    rxns_to_rule[rxn] = rx
    if rx == "":
        rxns_to_rule[rxn] = None
        continue
    expr = parse_expr(rx.replace("not", "~").replace("and", "&").replace("or", "|"))
    for x in expr.free_symbols:
        if str(x) in bnum_to_gene.keys():
            rx = rx.replace(str(x), bnum_to_gene[str(x)])
        else:
            rx = rx.replace(str(x), "True")    
    expr = parse_expr(rx)
    rxns_to_rule[rxn] = expr

print(rxns_to_rule)

#making a dictionary of reactions to vars and initial bounds
rxns_to_var_and_bounds = {}
for rxn in rxn_list:
    #get the reaction from the model
    rxn_obj = model.reactions.get_by_id(rxn)
    #get the bounds
    lb = rxn_obj.lower_bound
    ub = rxn_obj.upper_bound
    if lb == 0 and ub == 0:
        rxns_to_var_and_bounds[rxn] = [0, 0, 0]
    else:
        rxns_to_var_and_bounds[rxn] = [1, lb, ub]


# class state:
#     def __init__(self, model, initial_substrate_concs, supply, initial_X, tim,objective):
#         self.model = model
#         self.supply = supply #dictionary from uptake metabolites to supply values
#         self.substrate_concs = initial_substrate_concs #dictionary from metabolite to a list of concentrations
#         self.X = []
#         self.X.append(initial_X)
#         self.del_T = tim # in mins
#         self.counter = 0
#         if objective is None:
#             model.objective = 'Biomass_Ecoli_core_w_GAM'
#         else:
#             model.objective = objective
#         self.genes_to_rules = 
            
        
#     def convert_S2f(self,S_c,X,del_T):
#         #S_c is the concentration of the substrate
#         #X is the biomass concentration
#         #del_T is the time interval
#         #returns the flux of the substrate
#         return (S_c/(X*del_T))
        
#     def update_model(self):
#         print("_"*50)
#         print("Time: ",self.counter*self.del_T)
#         print("_"*50)
#         for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
#             S_c = conc[self.counter]
#             flux = self.convert_S2f(S_c,self.X[self.counter],self.del_T)
#             if model.reactions.get_by_id(rxn_id).upper_bound < -flux:
#                 model.reactions.get_by_id(rxn_id).upper_bound = 0
#             else:
#                 model.reactions.get_by_id(rxn_id).lower_bound = -flux
#                 model.reactions.get_by_id(rxn_id).upper_bound = 0
#         soln = model.optimize()
#         #extracting the values now
#         mu = soln.objective_value
#         d1 = {}
#         for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
#             d1[rxn_id] = soln.fluxes[rxn_id]
        
#         new_X = self.X[self.counter]*(np.exp(mu*(self.del_T/60)))
#         for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
#             self.substrate_concs[rxn_id].append(self.substrate_concs[rxn_id][self.counter] + (-d1[rxn_id]/mu)*(self.X[self.counter])*(1-np.exp(mu*(self.del_T/60))))
#         self.X.append(new_X)      
#         self.counter+=1
        
    
# # Define the initial conditions for the substrates and biomass
# glc_rx = 'EX_glc_e_'
# initial_substrate_concs = {glc_rx: [10]}  # example concentrations
# initial_biomass = 0.03  # example initial biomass
# timesteps = 300
# # print(model.objective)
# # soln1 = model.optimize()
# # print(soln1)
# model_state = state(model, initial_substrate_concs, None, initial_biomass, 1, "Biomass_Ecoli_core_N__w_GAM_")
# for i in range(timesteps):
#     model_state.update_model()
    
# X = np.arange(0,timesteps+1,1)
# Y = model_state.substrate_concs[glc_rx]
# print(Y)
# # plt.plot(X,Y)
# Y1 = model_state.X
# # plt.plot(X,Y1)
# plt.plot(Y1,Y)
# plt.show()
