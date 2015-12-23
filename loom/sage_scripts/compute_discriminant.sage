import sys
import pdb

#K.<I> = QuadraticField(-1)
#R.<z> = K[]
#S.<x> = R[]

#f_str = sys.argv[1]
#f = eval(f_str)
#delta = factor(f.polynomial(x).discriminant())

R.<z> = QQ[]
S.<x> = R[]

f_str = sys.argv[1]
f = eval(f_str)
delta = f.polynomial(x).discriminant()

print delta
