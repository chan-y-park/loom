import sys
from pprint import pprint

def translate_weight_data(data):
    """
    Weights are passed as instances of a certain sage class
    then just calling 'eval' on the answer will produce a 
    wrong result.
    This is because Sage deals with fractions while Python 
    does not.
    For example a weight [1/2, 1/2, -1/2] would be evaluated to
    [0,0,-1] wrongly.
    To fix this, we pick the elements one by one and approximate 
    them with reals.
    """

    ### Weight data is given as a dictionary with the following format
    ### {... , weight : multiplicity , ...}
    keys = data.keys()
    mults = data.values()

    weights = []
    for i in range(len(data) - 1):
        weight = map(float, list(keys[i].to_vector()))
        ### Now add the new entry
        weights.append(weight)

    last_weight =  map(float, list(keys[-1].to_vector()))
    ### Now add the last entry
    weights.append(last_weight)

    return [weights, mults]

algebra = sys.argv[1]
highest_weight_str = sys.argv[2]
highest_weight_list = eval(highest_weight_str)

ct = CartanType(algebra)

wcr = WeylCharacterRing(algebra)

fund_wts = wcr.fundamental_weights()

### Convert the highest weight from the base of coroots to the base of the ambient space.
highest_weight_internal = sum([x * fund_wts[i+1] for i, x in enumerate(highest_weight_list)])

rep = wcr(highest_weight_internal)

weights_mults = translate_weight_data(rep.weight_multiplicities())

pprint(weights_mults)