from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr
import json
# Define the symbols
# IclR, ArcA, FruR = symbols('IclR ArcA FruR')

# Define the string
s = "b1854  and  b1676"
a = s
#Replace English words with Python equivalents
s = s.replace("not", "~").replace("and", "&").replace("or", "|")
gene_to_bnum = {}
bnum_to_gene = {}
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/gene_to_bnum.json', 'r') as fp:
    gene_to_bnum = json.load(fp)
with open('/home/yashjonjale/Documents/iGEM/iGEM-IITB/FBA/rFBA/models/bnum_to_gene.json', 'r') as fp:
    bnum_to_gene = json.load(fp)
#Parse the string
expr = parse_expr(s)

print(bnum_to_gene)


print(expr)
print(expr.free_symbols)
for x in expr.free_symbols:
    if str(x) in bnum_to_gene.keys():
        s = s.replace(str(x), bnum_to_gene[str(x)])
    else:
        s = s.replace(str(x), "True")


print(s)
print(a)

d = "True  and  b1676"
expr = parse_expr(d)

print(expr)