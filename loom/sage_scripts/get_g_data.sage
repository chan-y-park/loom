import sys
import numpy

from pprint import pprint


def pick_basis(ffr_weights, algebra_type, algebra_rank):
    vector_list = [v.to_vector() for v in ffr_weights]
    basis = []

    if algebra_type == 'A':
        ### All the weights are linearly independent
        ### in the orthonormal basis.
        return vector_list

    elif algebra_type == 'D':
        ### Return a set of linearly independent weights
        ### whose coordinates in the orthonormal basis
        ### are positive.
        return vector_list[:algebra_rank]

    else:
        dim = algebra_rank
        
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
            #return map(list, basis)
            return basis
        
        else:
            raise RuntimeError(
                'pick_basis() '
                'This is not a complete basis.'
            )


def argsort(seq, reverse=False):
    #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    #by unutbu
    return sorted(range(len(seq)), key=seq.__getitem__, reverse=reverse)


def sort_ffr_weights(ffr_weights, algebra_type, algebra_rank):
    vector_list = [v.to_vector() for v in ffr_weights]
    if algebra_type == 'A':
        """
        Order the 1st fundamental weights in the following form,
        [(1, 0, 0, 0, 0, 0),
         (0, 1, 0, 0, 0, 0),
         ...,
         (0, 0, 0, 0, 0, 1)].
        """
        sorted_ffr_weights = [
            ffr_weights[i] for i in argsort(vector_list, reverse=True)
        ]
        return sorted_ffr_weights
    elif algebra_type == 'D':
        """
        Order the 1st fundamental weights in the following form,
        [(1, 0, 0, 0, 0),
         (0, 1, 0, 0, 0),
         ...,
         (0, 0, 0, 0, 1),
         (-1, 0, 0, 0, 0),
         ...
         (0, 0, 0, 0, -1)]      
        """
        sorted_ffr_positive_weights = [
            ffr_weights[i] 
            for i in argsort(vector_list, reverse=True)[:algebra_rank]
        ]
        sorted_ffr_negative_weights = [-v for v in sorted_ffr_positive_weights]
        return sorted_ffr_positive_weights + sorted_ffr_negative_weights
    elif algebra_type == 'E':
        return ffr_weights
    else:
        print('sort_ffr_weights(): there is no sorting rule for {}-type, '
              'returns without sorting. '.format(algebra_type))
        return ffr_weights


def main(root_system, highest_weight):
    algebra_type = root_system[0]
    algebra_rank = int(root_system[1:])
    R = RootSystem(root_system)
    ### A represent the weight space in the orthonormal basis.
    A = R.ambient_space()

    if algebra_type == 'A' or algebra_type == 'D':
        ### omega_1 is the 1st fundamental weight
        if root_system == 'D2':
            omega_1 = A.fundamental_weight(1) + A.fundamental_weight(2)
        else:
            omega_1 = A.fundamental_weight(1)   

    elif algebra_type == 'E':
        if algebra_rank == 6:
            ### omega_1 is the 1st fundamental weight
            omega_1 = A.fundamental_weight(1)   
        elif algebra_rank == 7:
            ### omega_1 is the 7th fundamental weight
            omega_1 = A.fundamental_weight(7)   
    else:
        raise RuntimeError(
            'get_g_data.sage.main(): '
            'There is no minuscule representation for this algebra.'
        )
    
    ### weyl_orbit_1 consists of weights in the 1st fundamental rep.
    weyl_orbit_1 = omega_1.orbit()
    ffr_weights = sort_ffr_weights(weyl_orbit_1, algebra_type, algebra_rank)

    #representation_index = eval(sys.argv[2])
    #if representation_index is not None:
    #    omega_n = A.fundamental_weight(representation_index)
    #    weyl_orbit_n = omega_n.orbit()
    #else:
    #    omega_n = None
    #    weyl_orbit_n = None

    # highet_weight is expressed in Dynkin labels.
    L = R.weight_space()
    L_fundamental_weights = L.fundamental_weights()
    L_highest_weight = 0
    for i, lambda_i in enumerate(highest_weight, start=1):
        L_highest_weight += lambda_i * L_fundamental_weights[i]
    if L_highest_weight == L.fundamental_weight(1):
        weights = ffr_weights
        multiplicities = [1] * len(ffr_weights)
    else:
        A_highest_weight = L_highest_weight.to_ambient()
        wcr = WeylCharacterRing(root_system)
        rep = wcr(A_highest_weight)
        ### Weight data is given as a dictionary with the following format
        ### {... , weight : multiplicity , ...}
        weight_multiplicities = rep.weight_multiplicities()
        weights = weight_multiplicities.keys()
        multiplicities = weight_multiplicities.values()
        
    ### Find a basis from the weights in the 1st fundamental rep,
    ### Then express the weights in the given rep in the basis.
    weight_basis = pick_basis(ffr_weights, algebra_type, algebra_rank)
    weight_coefficients = [Matrix(weight_basis).solve_left(weight.to_vector())
                           for weight in weights]

    data = {
        'ffr_weights': ffr_weights,
        'roots': A.roots(),
        'positive_roots': A.positive_roots(),
        'weights': weights,
        'multiplicities': multiplicities,
        'weight_basis': weight_basis,
        'weight_coefficients': weight_coefficients,
    }

    return data


root_system = sys.argv[1]
highest_weight = eval(sys.argv[2])

pprint(main(root_system, highest_weight))
