# import pdb


def get_A_D_curve_string(casimir_differentials, N):
    # g_type == 'A' or g_type == 'D':
    curve_str = 'x^{} '.format(N)

    for k, phi_k in casimir_differentials.items():
        curve_str += ' + ({}) '.format(phi_k)
        if k != N:
            curve_str += ' * x^{}'.format(N - k)

    return curve_str


# The E_6 curve in Class S form
# NOTE: here x really stands for \lambda
# NOTE: this being said, it's just a relabeling because \lambda = x dz/z
# this is terribly cofusing notation, should switch to writing
# curves entirely in terms of \lambda
# NOTE: The above comment applies only for pure gauge theories,
# where \lambda = v dt/t and therefore x = v/t.
# In general we can say \lambda = x dz, so x and \lambda can be used
# interchangeably without confusion when factoring out every term
# by dz^N, which is a common factor.
def get_E_6_curve_string(casimir_differentials):
    phi_dict = {}
    for k, phi_k in casimir_differentials.items():
        phi_dict['phi_{}'.format(k)] = phi_k

    p_1 = (
        '78*x^10 + 60*({phi_2})*x^8 + 14*({phi_2})^2*x^6'
        ' - 33*({phi_5})*x^5'
        ' + 2*({phi_6})*x^4 - 5*({phi_2})*({phi_5})*x^3'
        ' - ({phi_8})*x^2 - ({phi_9})*x - ({phi_5})^2'
    ).format(**phi_dict)
    p_2 = (
        '12*x^10 + 12*({phi_2})*x^8 + 4*({phi_2})^2*x^6'
        ' - 12*({phi_5})*x^5 + ({phi_6})*x^4'
        ' - 4*({phi_2})*({phi_5})*x^3 - 2*({phi_8})*x^2'
        ' + 4*({phi_9})*x + ({phi_5})^2'
    ).format(**phi_dict)
    q_1 = (
        '270*x^(15) + 342*({phi_2})*x^(13) + 162*({phi_2})^2*x^(11)'
        ' - 252*({phi_5})*x^(10) + (26*({phi_2})^3 + 18*({phi_6}))*x^9'
        ' - 162*({phi_2})*({phi_5})*x^8 + (6*({phi_2})*({phi_6})'
        ' - 27*({phi_8}))*x^7'
        ' - (30*({phi_2})^2*({phi_5}) - 36*({phi_9}))*x^6'
        ' + (27*({phi_5})^2 - 9*({phi_2})*({phi_8}))*x^5'
        ' - (3*({phi_5})*({phi_6}) - 6*({phi_2})*({phi_9}))*x^4'
        ' - 3*({phi_2})*({phi_5})^2*x^3 - 3*({phi_5})*({phi_9})*x'
        ' - ({phi_5})^3'
    ).format(**phi_dict)

    q_2 = (
        '1/(2*x^3)*(({q_1})^2 - ({p_1})^2*({p_2}))'
    ).format(**{'p_1': p_1, 'p_2': p_2, 'q_1': q_1})

    curve_string = (
        '(1/2)*(x^3)*({phi_12})^2 - ({q_1})*({phi_12}) + ({q_2})'
    ).format(**{'q_1': q_1, 'q_2': q_2, 'phi_12': phi_dict['phi_12']})

    return curve_string


def get_E_7_curve_string(casimir_differentials):
    phi_dict = {}
    for k, phi_k in casimir_differentials.items():
        phi_dict['phi_{}'.format(k)] = phi_k

    p_1 = (
        '1596*x^(10) + 88*({phi_2})^2*x^6 + 7*({phi_6})*x^4 '
        '+ 660*({phi_2})*x^8 - 2*({phi_8})*x^2 - ({phi_10})'
    ).format(**phi_dict)

    p_2 = (
        '16872*x^(15) + 11368*({phi_2})*x^(13) + 2952*({phi_2})^2*x^(11) '
        '+ (176*({phi_6}) + 264*({phi_2})^3)*x^9 '
        '+ (-100*({phi_8}) + (100/3)*({phi_2})* ({phi_6}))*x^7 '
        '+ (-(68/3)*({phi_2})*({phi_8}) + 68*({phi_10}))*x^5 '
        '+ ((2/9)*({phi_6})^2 -(4/3)*({phi_12}))*x^3 '
        '+ ((218/8613)*({phi_2})*({phi_6})^2 - (4/3)*({phi_14}))*x'
    ).format(**phi_dict)

    p_3 = (
        '44560*x^(20) + 41568*({phi_2})*x^(18) + 16080*({phi_2})^2*x^(16) '
        '+ (2880*({phi_2})^3 + (2216/3)*({phi_6}))*x^(14) '
        '+ (312*({phi_2})*({phi_6}) + 192*({phi_2})^4 - (1552/3)*({phi_8}))'
        '   *x^(12) '
        '+ (32*({phi_2})^2*({phi_6}) - 40*({phi_2})*({phi_10}) - '
        '   (64/3)*({phi_12}) + (11/3)*({phi_6})^2)*x^8 '
        '+ (-(416/3)*({phi_14}) - 16*({phi_2})^2*({phi_10}) '
        '   - (4/9)*({phi_6})*({phi_8}) - (32/9)*({phi_2})*({phi_12}) '
        '   + (27776/8613)*({phi_2})*({phi_6})^2)*x^6 '
        '+ ((3488/8613)*({phi_2})^2*({phi_6})^2 + (4/9)*({phi_8})^2 '
        '   - (64/3)*({phi_2})*({phi_14}) - (2/3)*({phi_6})*({phi_10}))*x^4 '
        '+ (4/3)*({phi_8})*({phi_10})*x^2 + ({phi_10})^2'
    ).format(**phi_dict)

    q = (
        '-28*x^(10) - (44/3)*({phi_2})*x^8 - (8/3)*({phi_2})^2*x^6 '
        '- (1/3)*({phi_6})*x^4 + (2/9)*({phi_8})*x^2 - (1/3)*({phi_10})'
    ).format(**phi_dict)

    r = (
        '148*x^(15) + 116*({phi_2})*x^(13) + 36*({phi_2})^2*x^(11) '
        '+ (4*({phi_2})^3 + (8/3)*({phi_6}))*x^9 '
        '+ ((2/3)*({phi_2})*({phi_6}) - 2*({phi_8}))*x^7 '
        '+ (2*({phi_10}) - (2/3)*({phi_2})*({phi_8}))*x^5 '
        '+ ((1/81)*({phi_6})^2 - (2/27)*({phi_12}))*x^3 '
        '+ ((109/8613)*({phi_2})*({phi_6})^2 - (2/3)*({phi_14}))*x'
    ).format(**phi_dict)

    sub_exprs = {'p_1': p_1, 'p_2': p_2, 'p_3': p_3, 'q': q, 'r': r}

    A_2 = '(9/16/x^2) * (6*({q})*({p_1}) - 3*({p_3}))'.format(**sub_exprs)

    A_1 = (
        '(9/16/x^2)^2 '
        '* (9*({q})^2*({p_1})^2 - 6*({r})*({p_1})*({p_2}) '
        '   - 12*({q})*({p_1})*({p_3}) + 3*({q})*({p_2})^2 + 3*({p_3})^2)'
    ).format(**sub_exprs)

    A_0 = (
        '-(9/16/x^2)^3 '
        '* (4*({r})^2*({p_1})^3 + 6*({q})*({r})*({p_1})^2*({p_2}) '
        '   + 9*({q})^2*({p_1})^2*({p_3}) - 6*({r})*({p_1})*({p_2})*({p_3}) '
        '   - 6*({q})*({p_1})*({p_3})^2 + 2*({r})*({p_2})^3 '
        '   + 3*({q})* ({p_2})^2*({p_3}) + ({p_3})^3)'
    ).format(**sub_exprs)

    curve_string = (
        '-(1/3^6)*x^2*(({phi_18})^3 + ({A_2})*({phi_18})^2 '
        '+ ({A_1})*({phi_18}) + ({A_0}))'
    ).format(
        **{'A_0': A_0, 'A_1': A_1, 'A_2': A_2, 'phi_18': phi_dict['phi_18']}
    )

    return curve_string


# TODO: make checks for the casimirs given by the user
# For example, I was able to insert the 0th Casimir with no complaint
# Moreover I didn get the curve I wanted because below we add +1 to phi_0
# We should prevent the wrong casimirs from being inserted at all.
def check_casimir_differentials(casimir_differentials, g_type, g_rank):
    if g_type == 'D':
        for k, phi_k in casimir_differentials.items():
            if k % 2 != 0:
                raise RuntimeError(
                    'check_casimir_differentials(): '
                    'phi_{} = {} dz^k '
                    'is an incorrect casimir differential '
                    'for a D-type curve.'
                    .format(k, phi_k)
                )


def get_ffr_curve_string(casimir_differentials, g_type, g_rank):
    """
    Construct a Seiberg-Witten curve in the 1st fundamental representation
    using Casimir differentials.
    """
    check_casimir_differentials(casimir_differentials, g_type, g_rank)

    if g_type == 'A':
        return get_A_D_curve_string(casimir_differentials, g_rank + 1)

    elif g_type == 'D':
        return get_A_D_curve_string(casimir_differentials, 2 * g_rank)

    elif g_type == 'E':
        if g_rank == 6:
            return get_E_6_curve_string(casimir_differentials)
        elif g_rank == 7:
            return get_E_7_curve_string(casimir_differentials)

    else:
        raise NotImplemented(
            'get_ffr_curve_string(): construction of a Seiberg-Witten curve '
            'of g = {}_{} is not implemented.'.format(g_type, g_rank)
        )
