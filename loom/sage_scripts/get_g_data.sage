import sys
from pprint import pprint
def pick_basis(vector_list)
    #vector_list = [vector(x) for x in eval(sys.argv[1])]
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
        return map(list, basis)
    
    else:
        raise ValueError('This is not a complete basis!')

root_system = sys.argv[1]

R = RootSystem(root_system)
A = R.ambient_space()

omega_1 = A.fundamental_weight(1)
weyl_orbit_1 = omega_1.orbit()

#representation_index = eval(sys.argv[2])
#if representation_index is not None:
#    omega_n = A.fundamental_weight(representation_index)
#    weyl_orbit_n = omega_n.orbit()
#else:
#    omega_n = None
#    weyl_orbit_n = None

# This is expressed in Dynkin labels.
highest_weight = eval(sys.argv[2])
L = R.weight_space()
L_fundamental_weights = L.fundamental_weights()
L_highest_weight = 0
for i, lambda_i in enumerate(highest_weight, start=1):
    L_highest_weight += lambda_i * L_fundamental_weights[i]
A_highest_weight = L_highest_weight.to_ambient()
wcr = WeylCharacterRing(root_system)
rep = wcr(A_highest_weight)
### Weight data is given as a dictionary with the following format
### {... , weight : multiplicity , ...}
weight_multiplicities = rep.weight_multiplicities()
weights = weight_multiplicities.keys(),
multiplicities = weight_multiplicities.values()

data = {
    #'omega_1': omega_1,
    #'omega_n': omega_n,
    'weyl_orbit_1': weyl_orbit_1,
    #'weyl_orbit_n': weyl_orbit_n,
    'roots': A.roots(),
    'positive_roots': A.positive_roots(),
    #'weight_multiplicities': weight_multiplicities,
    'weights': weights,
    'multiplicities': multiplicities,
}

pprint(data)
