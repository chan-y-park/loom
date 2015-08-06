#import subprocess
#import ast


#def positive_roots(algebra):
#    """
#    Arguments must be formatted as follows:
#    algebra = 'D3'
#    """
#    out = subprocess.check_output(
#        ["sage", "./sage_scripts/positive_roots.sage", 
#         algebra]
#    )
#    positive_roots = map(list, eval(out))
#    return positive_roots



#def weight_coefficients(weight, basis):
#    """
#    Expands weight in the elements of the basis provided.
#    The basis must not be overcomplete.
#    The format for weight is a list e.g. [1,0,0].
#    Similarly for basis [..., [0,1,-1], ...].
#    """
#    
#    out = subprocess.check_output(
#        ["sage", "./sage_scripts/weight_coefficients.sage", 
#         str(weight), str(basis)]
#    )
#    coefficients = eval(out)
#    return coefficients


#def pick_basis(vector_list):
#    """
#    Chooses a complete basis from the list of vectors.
#    """
#    out = subprocess.check_output(
#        ["sage", "./sage_scripts/pick_basis.sage", 
#         str(vector_list)]
#    )
#    basis = eval(out)
#    return basis
