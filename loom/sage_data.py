import subprocess
import ast

def weight_system(algebra, highest_weight):
    """
    Arguments must be formatted as follows:
    algebra = 'D3'
    highest_weight = [0, 1, 1]
    """
    out = subprocess.check_output(
        ["sage", "./sage_scripts/weight_system.sage", 
         algebra, 
         str(highest_weight)
        ]
    )
    weights, multiplicities = eval(out)
    return weights, multiplicities


def positive_roots(algebra):
    """
    Arguments must be formatted as follows:
    algebra = 'D3'
    """
    out = subprocess.check_output(
        ["sage", "./sage_scripts/positive_roots.sage", 
         algebra]
    )
    positive_roots = map(list, eval(out))
    return positive_roots

