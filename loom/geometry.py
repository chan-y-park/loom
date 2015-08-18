import sympy
import numpy
import logging
import pdb
from warnings import warn
from pprint import pformat

import sage_subprocess
from misc import ctor2, r2toc, get_root_multiplicity, PSL2C, n_nearest_indices

x, z = sympy.symbols('x z')


class GData:
    """
    ffr_weights:
        List of weights of the FIRST fundamental representation,
            [v_0, ... , v_k , ...],
        where 'v_k' are numpy arrays corresponding to weights.
        - For g=A_n Lie algebras, the weights are given in IR^{n+1} as
            v_0 = (1,0,...,0),
            v_1 = (0,1,0,..), 
            ...,
            v_n = (0,...,0,1). 
          In this case, it does not matter how we identify weights
          with sheets, since the Weyl group acts by permuting all of 
          them freely.
        - For g=D_n, the weights are given in IR^{n} as
            v_0 = (1,0,...,0),
            v_1 = (0,1,...,0),
            v_{n-1} = (0,...,0,1),
            v_n = (-1,0,...,0),
            v_{n+1} = (0,-1,...,0),
            v_{2n-1} = (0,...,0,-1).
          In this case, we diivde the sheets into positive and negative ones,
          and assign the weights accordingly.
          The assignment of positive sheets is almost arbitrary: from each pair
          of positive/negative sheets one can pick either, as long as one makes
          an even number of "sign mistakes". We don't keep track of this,
          as a result there is an ambiguity in distinguishing one spinor 
          representation from the other
        - For other algebras, the order is imposed by SAGE.
    """
    def __init__(self, root_system=None, representation_str=None):
        self.root_system = root_system
        ### type is 'A', 'D', or 'E'.
        self.type = root_system[0]
        self.rank = eval(root_system[1:])

        representation = eval(representation_str)
        if isinstance(representation, int):
            ### Representation is specified as an index
            ### of a fundamental representation, i.e. n of \omega_n.
            self.fundamental_representation_index = representation
            self.highest_weight = [
                1 if i == (self.fundamental_representation_index - 1) else 0
                for i in range(self.rank)
            ]
        elif isinstance(representation, list):
            ### Representation is specified in the coroot(Dynkin) basis.
            self.highest_weight = representation
            height = 0
            for i, n_i in enumerate(self.highest_weight):
                height += n_i
                if n_i == 1:
                    self.fundamental_representation_index = i + 1
            if height > 1:
                self.fundamental_representation_index = None 
                warn('{} is not a fundamental representation.'
                     .format(representation_str))
                
        sage_data = sage_subprocess.get_g_data(
            root_system, 
            self.highest_weight,
        )
        logging.info("Representation data of ({}, {}) retrieved from SAGE:\n"
                     "{}".format(root_system, representation_str,
                                 pformat(sage_data)))

        self.ffr_weights = numpy.array(sage_data['ffr_weights'])
        self.roots = numpy.array(sage_data['roots'])
        self.positive_roots = numpy.array(sage_data['positive_roots'])
        self.weights = numpy.array(sage_data['weights'])
        self.multiplicities = numpy.array(sage_data['multiplicities'])
        self.weight_basis = numpy.array(sage_data['weight_basis']),
        ### The i-th row of self.coefficients is the representation
        ### of self.weights[i] in the self.basis.
        self.weight_coefficients = numpy.array(sage_data['weight_coefficients'])

    def ordered_weight_pairs(self, root, ffr=False):
        pairs = []

        if ffr==False:
            for i, w_1 in enumerate(self.weights):
                for j, w_2 in enumerate(self.weights):
                    if numpy.array_equal(w_2 - w_1, root):
                        pairs.append([i, j])
        elif ffr==True:
            for i, w_1 in enumerate(self.ffr_weights):
                for j, w_2 in enumerate(self.ffr_weights):
                    if numpy.array_equal(w_2 - w_1, root):
                        pairs.append([i, j])

        return pairs

    def weyl_monodromy(self, root, bp, direction):
        ### TODO ###
        ### For now, we ASSUME that branch points are 
        ### of SQUARE-ROOT type. Need to update this 
        ### for more general cases.
        bp_root = bp.positive_roots[0]

        new_root = root - (2 * bp_root * (numpy.dot(root, bp_root))
                           / numpy.dot(bp_root,bp_root))

        return new_root


class RamificationPoint:
    def __init__(self, z=None, x=None, i=None, label=None, is_puncture=False):
        self.z = z
        self.x = x
        self.i = i
        self.label = label
        self.is_puncture = is_puncture

    def __str__(self):
        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'x': ctor2(self.x),
            'i': ctor2(self.i),
            'label': self.label,
            'is_puncture': self.is_puncture
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.x = r2toc(json_data['x'])
        self.i = r2toc(json_data['i'])
        self.label = json_data['label']
        self.is_puncture = json_data['is_puncture']


class PuncturePoint:
    def __init__(self, z, cutoff, label):
        self.z = z
        self.cutoff = cutoff
        self.label = label

    def __eq__(self, other):
        return self.label == other.label


class SWCurve:
    def __init__(self, differentials=None, Cz=None, parameters=None,
                 ffr_weights=None,):
        if ffr_weights is not None:
            ### Build a cover in the first fundamental representation.
            N = len(ffr_weights)
            ffr_eq_str = 'x^{} '.format(N)
            for k, u_k in differentials.iteritems():
                ffr_eq_str += '+ ({}) '.format(u_k)
                if k != N:
                    ffr_eq_str += '* x^{}'.format(N-k)

            self.sym_eq = sympy.simplify(
                sympy.sympify(ffr_eq_str).subs(z, Cz)
            )
            self.num_eq = self.sym_eq.subs(parameters)
        else:
            ### TODO: Need to build a cover in a general representation
            ### from the differentials, using symmetric polynomials.
            logging.warning(
                'class SWCurve with a general representation '
                'is not implemented yet.'
            )
            self.sym_eq = None 
            self.num_eq = None


    def get_ramification_points(self, accuracy, punctures=None):
        f = self.num_eq
        if f is None:
            raise NotImplementedError

        ramification_points = []

        # NOTE: solve_poly_system vs. solve
        #sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
        #sols = sympy.solve([f, f.diff(x)], z, x)
        #if sols is None:
        #    # Use Sage instead
        #    sols = sage_subprocess.solve_poly_system([f, f.diff(x)])
        sols = sage_subprocess.solve_poly_system([f, f.diff(x)])
        for z_0, x_0 in sols:
            if (len(punctures) > 0 and
                (min([abs(z_0 - p) for p in punctures]) < accuracy)
            ):
                continue
            fx_at_z_0 = f.subs(z, z_0)
            fx_at_z_0_coeffs = map(
                complex, sympy.Poly(fx_at_z_0, x).all_coeffs()
            )
            mx = get_root_multiplicity(
                fx_at_z_0_coeffs, complex(x_0), accuracy
            )
            if mx > 1:
                fz_at_x_0 = f.subs(x, x_0)
                fz_at_x_0_coeffs = map(complex, 
                                       sympy.Poly(fz_at_x_0, z).all_coeffs())
                
                mz = get_root_multiplicity(fz_at_x_0_coeffs, complex(z_0),
                                           accuracy)
                if mz > 1:
                    continue
                label = ('ramification point #{}'
                         .format(len(ramification_points)))
                rp = RamificationPoint(complex(z_0), complex(x_0), mx, label)
                logging.info("{}: z = {}, x = {}, i = {}."
                             .format(label, rp.z, rp.x, rp.i))
                ramification_points.append(rp)
        return ramification_points


    def get_xs(self, z_0):
        """
        Return a numpy array of x-coordinates over z = z_0.
        """
        if self.num_eq is None:
            raise NotImplementedError

        fx = self.num_eq.subs(z, z_0)
        #xs = sympy.solve(fx, x)
        ### The following may fail when the curve is not a polynomial.
        sym_poly = sympy.Poly(fx, x, domain='CC')
        coeff_list = map(complex, sym_poly.all_coeffs())
        return numpy.roots(coeff_list)


class SWDiff:
    def __init__(self, v_str, g_data=None, Cz=None, dCz=None, parameters=None):
        ### sym_v is a SymPy expression. 
        self.sym_v = sympy.simplify(
            sympy.sympify(v_str).subs(z, Cz) * dCz
        )
        ### num_v is from sym_v with its parameters 
        ### substituted with numerical values.
        self.num_v = self.sym_v.subs(parameters)


class SWData(object):
    """
    A class containing a Seiberg-Witten curve in the first
    fundamental representation,
        \lambda^N + \sum_{k=2}^N \phi_k \lambda^{N-k}, 
    where \lambda is the Seiberg-Witten differential of the form 
        \lambda = x dz
    and \phi_k = u_k(z) dz^k, and another curve in the representation
    given in the configuration.
    """
    def __init__(self, config):
        self.punctures = None
        self.parameters = config['sw_parameters']
        self.differentials = eval(config['differentials'])
        self.g_data = GData(config['root_system'], config['representation'])

        ### PSL2C-transformed z & dz
        Cz = PSL2C(config['mt_params'], z, inverse=True) 
        dCz = Cz.diff(z)

        self.ffr_curve = SWCurve(
            differentials=self.differentials, 
            Cz=Cz,
            parameters=self.parameters,
            ffr_weights=self.g_data.ffr_weights,
        )
        logging.info('Seiberg-Witten curve in the 1st fundamental '
                     'representation: {} = 0'
                     .format(sympy.latex(self.ffr_curve.num_eq)))
        ### TODO: SWCurve in a general representation.
        if self.g_data.fundamental_representation_index == 1:
            self.curve = self.ffr_curve
        else:
            logging.warning(
                'Seiberg-Witten curve in a general representation '
                'is not implemented yet.'
            )


        ### Seiberg-Witten differential
        self.diff = SWDiff(
            'x',
            g_data=self.g_data,
            Cz=Cz,
            dCz=dCz,
            parameters=self.parameters,
        )

        logging.info('\nSeiberg-Witten differential: %s dz\n',
                     sympy.latex(self.diff.num_v))

        if config['punctures'] is None:
            self.punctures = []
        else:
            self.punctures = [
                PSL2C(config['mt_params'], p, numerical=True)
                for p in config['punctures']
            ]

        logging.info(
            "Calculating ramification points of the Seiberg-Witten curve "
            "in the first fundamental rep."
        )
        self.ffr_ramification_points = self.ffr_curve.get_ramification_points(
            config['accuracy'], punctures=self.punctures,
        )


    def get_aligned_xs(self, z_0):
        ### Returns (aligned_ffr_xs, aligned_xs), where each element is
        ### a numpy array of x-coordinates of the fibers over z.
        ### The order of x's is the same as the order of the weights
        ### in g_data.weights.

        algebra_type = self.g_data.type
        algebra_rank = self.g_data.rank
        fund_rep_index = self.g_data.fundamental_representation_index

        ### First order x's of the first fundamental cover
        ### according to the order of ''g_data.weights'',
        ### then construct the list of x's for the given representation
        ### ordered according to weights.
        ffr_xs = self.ffr_curve.get_xs(z_0) 

        if algebra_type == 'A':
            """
            Can consider ffr_xs to be aligned according to ffr_weights,
            [(1, 0, 0, 0, 0, 0),
             (0, 1, 0, 0, 0, 0),
             ...,
             (0, 0, 0, 0, 0, 1)].
            """
            aligned_ffr_xs = ffr_xs

            if fund_rep_index == 1:
                xs = aligned_ffr_xs
            else:
                xs = self.get_xs_of_weights_from_ffr_xs(aligned_ffr_xs)

        elif algebra_type == 'D':
            """
            Align ffr_xs according to ffr_weights,
            [(1, 0, 0, 0, 0),
             (0, 1, 0, 0, 0),
             ...,
             (0, 0, 0, 0, 1),
             (-1, 0, 0, 0, 0),
             ...
             (0, 0, 0, 0, -1)]      
            """
            sorted_ffr_xs = numpy.array(
                sorted(ffr_xs, key=lambda z: (z.real, z.imag), reverse=True,)
            )
            ### Pick x's corresponding to the positive weights.
            ### The order among the positive x's is arbitrary.
            positive_xs = sorted_ffr_xs[:algebra_rank]
            ### Then pick an x corresponding to each negative weight
            ### aligned according to the positive x's.
            negative_xs = numpy.zeros_like(positive_xs)
            for nx in sorted_ffr_xs[algebra_rank:]:
                difference = numpy.fromiter(
                    (abs(px - (-nx)) for px in positive_xs),
                    dtype=float,
                )
                j = difference.argsort()[0]
                ### Check the pairing of positive and negative x's.
                px_j = positive_xs[j]
                if numpy.isclose(px_j, -nx) is False:
                    warn("get_ordered_xs(): No pairing of x's in the D-type, "
                         "({}, {}) != (x, -x).".format(px_j, nx))
                else:
                    ### Put the negative x at the same index
                    ### as its positive pair.
                    negative_xs[j] = nx
            aligned_ffr_xs = numpy.concatenate((positive_xs, negative_xs))

            if fund_rep_index == 1:
                xs = aligned_ffr_xs
            else:
                xs = self.get_xs_of_weights_from_ffr_xs(aligned_ffr_xs)

        elif g_data.type == 'E':
            ### !!!!!
            ### Here I am pairing any sheet with any weight of E_6 and E_7 
            ### However, the Weyl group does not contains permutations 
            ### of 27 (resp 56) elements so it's probably NOT OK to 
            ### pair sheets with weights as we like.
            ### !!!!!
            #xs = ffr_xs
            raise NotImplementedError
     
        return (aligned_ffr_xs, xs)


    def get_xs_of_weights_from_ffr_xs(self, ffr_xs):
        g_data = self.g_data
        xs = numpy.zeros(len(g_data.weights), dtype=complex)

        if g_data.type == 'A' or g_data.type == 'D':
            for i, cs in enumerate(g_data.weight_coefficients):
                for j, c_j in enumerate(cs):
                    xs[i] += c_j * ffr_xs[j]
        else:
            raise NotImplementedError

        return xs


def find_xs_at_z_0(sw_data, z_0, x_0=None, num_x=1, ffr=False):
    """
    Get x's above z_0 and return the num_x of them 
    which are nearest to x_0.
    """
    
    xs_at_z_0 = sw_data.get_sheets_at_z(z_0, ffr=ffr).values()
    if x_0 is None:
        return xs_at_z_0
    else:
        return sorted(xs_at_z_0,
                      lambda x1, x2: cmp(abs(x1 - x_0), abs(x2 - x_0)))[:num_x]
   

def get_local_sw_diff(sw, ramification_point, ffr=None):
    if ffr==None or ffr==False:
        rp = ramification_point
        num_eq = sw.curve.num_eq
        num_v = sw.diff.num_v
    elif ffr==True:
        rp = ramification_point
        num_eq = sw.ffr_curve.num_eq
        num_v = sw.diff.num_v

    ### use Dz = z - rp.z & Dx = x - rp.x
    Dz, Dx = sympy.symbols('Dz, Dx')
    local_curve = (
        num_eq.subs(x, rp.x+Dx).subs(z, rp.z+Dz)
        .series(Dx, 0, rp.i+1).removeO()
        .series(Dz, 0, 2).removeO()
    )
    ### curve_at_rp = a(z - rp.z) + b(x - rp.x)^(rp.i)
    a = local_curve.n().coeff(Dz).coeff(Dx, 0)
    b = local_curve.n().coeff(Dx**rp.i).coeff(Dz, 0)
    ### Dx = Dx(Dz)
    Dx_Dz = (-(a/b)*Dz)**sympy.Rational(1, rp.i)
    local_diff = (
        num_v.subs(x, rp.x+Dx_Dz).subs(z, rp.z+Dz)
        .series(Dz, 0, 1).removeO()
    )
    # get the coefficient and the exponent of the leading term
    (diff_c, diff_e) = local_diff.leadterm(Dz)
    if diff_e == 0:
        # remove the constant term from the local_diff
        local_diff -= local_diff.subs(Dz, 0)
        (diff_c, diff_e) = local_diff.leadterm(Dz)

    return (complex(diff_c.n()), diff_e)
