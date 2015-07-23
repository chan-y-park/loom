import sys
from pprint import pprint

vector_list = [vector(x) for x in eval(sys.argv[1])]
dim = len(vector_list[0])

basis = []
r = Matrix(basis).rank()

for v in vector_list:
    new_basis = basis + [v]
    new_r = Matrix(new_basis).rank()

    if new_r > r:
        basis = new_basis
        r = new_r

    else:
        pass

if r == dim:    
    pprint(map(list, basis))

else:
    raise ValueError('This is not a complete basis!')