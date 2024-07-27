import numpy
import config

well_labels = []
for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
    for number in range(1,13):
        well_labels.append('%s%d' % (letter , number))


antibiotic_labels = ['cN',  'c11',  'c9',  'c7']
antibiotic_conc = 25*numpy.asarray([0, 2**-11, 2**-9,  2**-7])

antibiotic_labels_dict = {}
for a_idx, a in enumerate(antibiotic_labels):
    antibiotic_labels_dict[a] = antibiotic_conc[a_idx]

phage_labels = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']
#phage_labels = phage_labels[::-1]
phage_dilution = numpy.asarray([10**(-1*int(s[1])) for s in phage_labels])
phage_labels_dict = {}
for a_idx, a in enumerate(phage_labels):
    phage_labels_dict[a] = phage_dilution[a_idx]



def parse_metadata():

    label_well_dict = {}

    line_count = 0
    for line in open('%slabels.csv' % config.data_directory,  encoding='utf-8-sig'):

        line = line.strip().split(',')

        for item_idx, item in enumerate(line):
            #label_well_dict[well_labels[(12*line_count) + item_idx]] = item
            label_well_dict[well_labels[(12*line_count) + item_idx]] = item

        line_count += 1

    return label_well_dict




def parse_data():

    data_dict = {}
    data_dict['hours'] = []

    label_well_dict = parse_metadata()

    for key, value in label_well_dict.items():
        data_dict[value] = []

    plate_data = open('%skitp_qbio_turner_2024_07_24_abx_phage_conc_gradient.csv' % config.data_directory, 'r',  encoding='utf-8-sig')
    header = plate_data.readline()
    header = header.strip().split(',')
    for line in plate_data:
        
        line = line.strip().split(',')
        line = [float(s) for s in line]
        data_dict['hours'].append(line[0]/(60*60))

        for od_idx, od in enumerate(line[1:]):
            
            data_dict[label_well_dict[header[od_idx+1]]].append(od)

    plate_data.close()

    for key, value in data_dict.items():

        data_dict[key] = numpy.asarray(value)

    return data_dict




def fd_weights_all(x, x0=0, n=1):
    """
    Return finite difference weights for derivatives of all orders up to n.
    Parameters
    ----------
    x : vector, length m
        x-coordinates for grid points
    x0 : scalar
        location where approximations are to be accurate
    n : scalar integer
        highest derivative that we want to find weights for
    Returns
    -------
    weights :  array, shape n+1 x m
        contains coefficients for the j'th derivative in row j (0 <= j <= n)
    Notes
    -----
    The x values can be arbitrarily spaced but must be distinct and len(x) > n.
    The Fornberg algorithm is much more stable numerically than regular
    vandermonde systems for large values of n.
    See also
    --------
    fd_weights
    References
    ----------
    B. Fornberg (1998)
    "Calculation of weights_and_points in finite difference formulas",
    SIAM Review 40, pp. 685-691.
    http://www.scholarpedia.org/article/Finite_difference_method
    """
    m = len(x)
    _assert(n < m, 'len(x) must be larger than n')

    weights = numpy.zeros((m, n + 1))
    _fd_weights_all(weights, x, x0, n)
    return weights.T


# from numba import jit, float64, int64, int32, int8, void
# @jit(void(float64[:,:], float64[:], float64, int64))
def _fd_weights_all(weights, x, x0, n):
    m = len(x)
    c1, c4 = 1, x[0] - x0
    weights[0, 0] = 1
    for i in range(1, m):
        j = numpy.arange(0, min(i, n) + 1)
        c2, c5, c4 = 1, c4, x[i] - x0
        for v in range(i):
            c3 = x[i] - x[v]
            c2, c6, c7 = c2 * c3, j * weights[v, j - 1], weights[v, j]
            weights[v, j] = (c4 * c7 - c6) / c3
        weights[i, j] = c1 * (c6 - c5 * c7) / c2
        c1 = c2


def fd_weights(x, x0=0, n=1):
    """
    Return finite difference weights for the n'th derivative.
    Parameters
    ----------
    x : vector
        abscissas used for the evaluation for the derivative at x0.
    x0 : scalar
        location where approximations are to be accurate
    n : scalar integer
        order of derivative. Note for n=0 this can be used to evaluate the
        interpolating polynomial itself.
    Examples
    --------
    >>> import numpy as np
    >>> import numdifftools.fornberg as ndf
    >>> x = numpy.linspace(-1, 1, 5) * 1e-3
    >>> w = ndf.fd_weights(x, x0=0, n=1)
    >>> df = numpy.dot(w, numpy.exp(x))
    >>> numpy.allclose(df, 1)
    True
    See also
    --------
    fd_weights_all
    """
    return fd_weights_all(x, x0, n)[-1]



def _assert(cond, msg):
    if not cond:
        raise ValueError(msg)


def fd_derivative(fx, x, n=1, m=2):
    """
    Return the n'th derivative for all points using Finite Difference method.
    Parameters
    ----------
    fx : vector
        function values which are evaluated on x i.e. fx[i] = f(x[i])
    x : vector
        abscissas on which fx is evaluated.  The x values can be arbitrarily
        spaced but must be distinct and len(x) > n.
    n : scalar integer
        order of derivative.
    m : scalar integer
        defines the stencil size. The stencil size is of 2 * mm + 1
        points in the interior, and 2 * mm + 2 points for each of the 2 * mm
        boundary points where mm = n // 2 + m.
    fd_derivative evaluates an approximation for the n'th derivative of the
    vector function f(x) using the Fornberg finite difference method.
    Restrictions: 0 < n < len(x) and 2*mm+2 <= len(x)
    Examples
    --------
    >>> import numpy as np
    >>> import numdifftools.fornberg as ndf
    >>> x = numpy.linspace(-1, 1, 25)
    >>> fx = numpy.exp(x)
    >>> df = ndf.fd_derivative(fx, x, n=1)
    >>> numpy.allclose(df, fx)
    True
    See also
    --------
    fd_weights
    """
    num_x = len(x)
    _assert(n < num_x, 'len(x) must be larger than n')
    _assert(num_x == len(fx), 'len(x) must be equal len(fx)')

    du = numpy.zeros_like(fx)

    mm = n // 2 + m
    size = 2 * mm + 2  # stencil size at boundary
    # 2 * mm boundary points
    for i in range(mm):
        du[i] = numpy.dot(fd_weights(x[:size], x0=x[i], n=n), fx[:size])
        du[-i - 1] = numpy.dot(fd_weights(x[-size:], x0=x[-i - 1], n=n), fx[-size:])

    # interior points
    for i in range(mm, num_x - mm):
        du[i] = numpy.dot(fd_weights(x[i - mm:i + mm + 1], x0=x[i], n=n),
                       fx[i - mm:i + mm + 1])

    return du

