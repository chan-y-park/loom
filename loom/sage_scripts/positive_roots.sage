#import sys

#NOTE: Don't need this, use __future__.division.
#def translate_root_data(data):
#    """
#    Roots are passed as instances of a certain sage class
#    then just calling 'eval' on the answer will produce a 
#    wrong result.
#    This is because Sage deals with fractions while Python 
#    does not.
#    For example a root [1/2, 1/2, -1/2] would be evaluated to
#    [0,0,-1] wrongly.
#    To fix this, we pick the elements one by one and approximate 
#    them with reals.
#    """
#
#    return [map(float, list(root.to_vector())) for root in data]
#
#root_system = sys.argv[1]
#
#Phi = RootSystem(root_system)
#
#print Phi.ambient_space().positive_roots()
