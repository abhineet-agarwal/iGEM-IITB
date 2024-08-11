import cobra

# Load the model
model = cobra.io.read_sbml_model("/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/Metabolic_Reconstructions/MODEL1108160000_uri.xml")

# Print the model - debug
# print(model)

# Perform FBA
solution = model.optimize()

# The solution

print(solution)

