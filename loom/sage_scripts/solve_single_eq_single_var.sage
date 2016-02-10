import sys
import pdb

#z = var('z')
#eqs = []

#precision = eval(sys.argv[1])
#for eq_str in sys.argv[2:]:
#    eqs.append(eval(eq_str))


#sols = solve(eqs, z)

#sols_str = []
#messages = []
## the following is a dummy variable here, 
## but kept for convenience of compatibility
## when switching between two different methods
#sols_mul = []
#for sol in sols:
#    z_ans= sol
#
#    try:
#        z_n = z_ans.right().n(digits=precision)
#    except TypeError as e:
#        for msg in e.args:
#            messages.append('solve_single_eq_z.sage: {}'.format(e))
#        messages.append('z = {}'.format(z_ans))
#        continue
#
#    z_re = z_n.real()
#    z_im = z_n.imag()
#
#    sols_str.append(
#        (str(z_re), str(z_im))
#    )
#    sols_mul.append(1)
#
#print (sols_str, sols_mul, messages)

### The following gives better numerics.
### But for computing the discriminant, 
### the above method is actually often better.
### This is because we only use the discriminant 
### Method for highly degenerate curves, and then 
### it can be solved exactly.
### However, when we want to use the discriminant 
### method for more general curves, we need to use 
### this routine.

precision = eval(sys.argv[1])
eq_str = sys.argv[2]
var_name = sys.argv[3]

var(var_name)
R = PolynomialRing(CDF, ('%s' % var_name))
eq = R(eval(eq_str))
sols = eq.roots()

sols_str = []
sols_mul = []
messages = []

for sol in sols:
    z_ans = sol[0]
    z_mul = sol[1]

    try:
        z_n = z_ans
    except TypeError as e:
        for msg in e.args:
            messages.append('solve_single_eq.sage: {}'.format(e))
        messages.append(var + ' = {}'.format(z_ans))
        continue

    z_re = z_n.real()
    z_im = z_n.imag()

    sols_str.append(
        (str(z_re), str(z_im))
    )
    sols_mul.append(z_mul)
print (sols_str, sols_mul, messages)

