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
    To fix this, we picke the elements one by one and approximate 
    them with reals.
    """

    keys = [w for w in data]
    mults = str([int(data[k]) for k in keys])
    weights = '['
    for i in range(len(data) - 1):
        weight = keys[i]
        l = len(list(eval(str(weight))))
        weight_as_list = [N(weight[i]) for i in range(l)]
        ### Now add the new entry
        weights += str(weight_as_list) + ', '

    last_weight = keys[-1]
    l = len(list(eval(str(weight))))
    last_weight_as_list = [N(last_weight[i]) for i in range(l)]
    ### Now add the last entry
    weights += str(last_weight_as_list) + ']'

    return '[' + weights + ', ' + mults + ']'
    # return mults

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