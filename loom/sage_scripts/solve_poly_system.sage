import sys
import pdb

z, x = var('z x')
poly_system = []
for poly_str in sys.argv[1:]:
    poly_system.append(eval(poly_str))
sols = solve(poly_system, z, x)

sols_str = []
for sol in sols:
    z_ans, x_ans = sol
    sols_str.append([complex(z_ans.right()), complex(x_ans.right())])
print sols_str
