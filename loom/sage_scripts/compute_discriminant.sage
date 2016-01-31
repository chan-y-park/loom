import sys
import pdb

#I = CC(0,1)
#R.<z> = CC[]
#S.<x> = R[]

#f_str = sys.argv[1]
#f = eval(f_str)
#delta = f.polynomial(x).discriminant()
#print delta


R.<x,z,I>=QQ[]
f_str = sys.argv[1]
f = eval(f_str)
delta = f.discriminant(x).subs(I = CC(0,1))
print delta

