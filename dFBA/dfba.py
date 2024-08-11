import cobra
import numpy as np
from matplotlib import pyplot as plt

# Load your model
model = cobra.io.read_sbml_model('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/Metabolic_Reconstructions/ecoli_core.xml')


class state:
    def __init__(self, model, initial_substrate_concs, supply, initial_X, tim,objective):
        self.model = model
        self.supply = supply #dictionary from uptake metabolites to supply values
        self.substrate_concs = initial_substrate_concs #dictionary from metabolite to a list of concentrations
        self.X = []
        self.X.append(initial_X)
        self.del_T = tim # in mins
        self.counter = 0
        if objective is None:
            model.objective = 'Biomass_Ecoli_core_w_GAM'
        else:
            model.objective = objective
            
        
    def convert_S2f(self,S_c,X,del_T):
        #S_c is the concentration of the substrate
        #X is the biomass concentration
        #del_T is the time interval
        #returns the flux of the substrate
        return (S_c/(X*del_T))
        
    def update_model(self):
        print("_"*50)
        print("Time: ",self.counter*self.del_T)
        print("_"*50)
        for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
            S_c = conc[self.counter]
            flux = self.convert_S2f(S_c,self.X[self.counter],self.del_T)
            if model.reactions.get_by_id(rxn_id).upper_bound < -flux:
                model.reactions.get_by_id(rxn_id).upper_bound = 0
            else:
                model.reactions.get_by_id(rxn_id).lower_bound = -flux
                model.reactions.get_by_id(rxn_id).upper_bound = 0
        soln = model.optimize()
        #extracting the values now
        mu = soln.objective_value
        d1 = {}
        for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
            d1[rxn_id] = soln.fluxes[rxn_id]
        
        new_X = self.X[self.counter]*(np.exp(mu*(self.del_T/60)))
        for i, (rxn_id, conc) in enumerate(self.substrate_concs.items()):
            self.substrate_concs[rxn_id].append(self.substrate_concs[rxn_id][self.counter] + (-d1[rxn_id]/mu)*(self.X[self.counter])*(1-np.exp(mu*(self.del_T/60))))
        self.X.append(new_X)      
        self.counter+=1
        
    
if __name__ == '__main__':
    # Define the initial conditions for the substrates and biomass
    glc_rx = 'EX_glc_e_'
    initial_substrate_concs = {glc_rx: [10]}  # example concentrations
    initial_biomass = 0.03  # example initial biomass
    timesteps = 300
    # print(model.objective)
    # soln1 = model.optimize()
    # print(soln1)
    model_state = state(model, initial_substrate_concs, None, initial_biomass, 1, "Biomass_Ecoli_core_N__w_GAM_")
    for i in range(timesteps):
        model_state.update_model()
        
    X = np.arange(0,timesteps+1,1)
    Y = model_state.substrate_concs[glc_rx]
    print(Y)
    # plt.plot(X,Y)
    Y1 = model_state.X
    # plt.plot(X,Y1)
    plt.plot(Y1,Y)
    plt.show()

