import sys
import pdb

R.<x> = PolynomialRing(CDF)
eqs = []

precision = eval(sys.argv[1])
for eq_str in sys.argv[2:]:
    eqs.append(R(eval(eq_str)))


#sols, mults = solve(eqs, x, multiplicities=True)
sols = eqs[0].roots()

sols_str = []
sols_mul = []
messages = []

for sol in sols:
    x_ans = sol[0]
    x_mul = sol[1]

    try:
        x_n = x_ans
    except TypeError as e:
        for msg in e.args:
            messages.append('solve_single_eq.sage: {}'.format(e))
        messages.append('x = {}'.format(x_ans))
        continue

    x_re = x_n.real()
    x_im = x_n.imag()

    sols_str.append(
        (str(x_re), str(x_im))
    )
    sols_mul.append(x_mul)


print (sols_str, sols_mul, messages)
#print (sols_str, messages)
