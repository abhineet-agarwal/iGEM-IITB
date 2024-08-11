import cobra
import scipy.io

mat_data = scipy.io.loadmat('your_matlab_model.mat')
model = cobra.io.load_matlab_model(mat_data, variable_name='model') 


solution = model.optimise()


# # Load the model
# # model = cobra.io.read_sbml_model("/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/Extended_FBA_ehux/S1_iEH410_model_files/iEH410.xml")
# model = cobra.io.load_matlab_model("/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/Extended_FBA_ehux/S1_iEH410_model_files/iEH410.mat",variable_name="model")

# solution = model.optimize()

print(solution.objective_value)
print("fluxes:  ")
print(solution.fluxes)