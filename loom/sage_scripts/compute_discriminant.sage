import sys
import pdb

R.<x,z> = QQ[]

f_str = sys.argv[1]
f = eval(f_str)
delta = factor(f.polynomial(x).discriminant())

print delta
