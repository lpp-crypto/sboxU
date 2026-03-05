"""It is possible to investigate to define and to compute the probability that a given S-box is "at least as good" or "at least as bad" as a random permutation or function, as explained in [AC:BonPerTia19].

This reasoning is based on the probabilities for the DDT and LAT given in [JMC:DaeRij07], while the probabilities from the BCT are from [AC:BonPerTia19].

"""

from sage.all import RealField, RealNumber, imag_part, exp, factorial, binomial, Infinity, Integer, floor
import itertools

from sboxUv2.core import Sb, is_permutation
from sboxUv2.statistics.cython_functions import differential_spectrum, walsh_spectrum, boomerang_spectrum


from sboxUv2.config import DEFAULT_HIGH_PRECISION



# !SECTION! Probability distributions


# !SUBSECTION! Functions handling probability distributions


def expected_maximum_distribution(
        in_length,
        out_length,
        proba_func,
        n_coeffs,
        vmin=2,
        vmax=100):
    """Assuming that a table contains 2**in_length-1 non trivial rows of 2**out_length-1 non-trivial columns, returns the probability distribution of the maximum coefficient in the table as dictionary.
    
    Args:
        in_length (int): the number of bits in the input
        out_length (int): the number of bits in the output
        proba_func (function): a function taking as input an input length, an output length and a coefficient, and which returns the probability that this coefficient occurs given the other parameters
        n_coeffs (int): the total number of non-trivial coefficients
        vmin (int): the minimum value of the differential uniformity to consider
        vmax (int): the maximum value of the differential uniformity to consider

    Returns:
        dict: A dictionary `d` such that `d[i]` is the probability that the maximum coefficient is equal to `i`.
    
    """
    cumul = {}
    m, n = in_length, out_length
    power = RealNumber(n_coeffs)
    for k in range(vmin, vmax+2):
        to_add = sum(proba_func(m, n, z) for z in range(0, k+1))
        to_add = to_add**power
        if to_add > 0:
            cumul[k] = to_add
        elif k-1 in cumul.keys():
            cumul[k] = cumul[k-1]
    result = {}
    for i, k in enumerate(sorted(cumul.keys())):
        if i > 0:
            result[k] = cumul[k] - cumul[k-1]
    return result


# !SUBSECTION! DDT

def ddt_coeff_probability(
        in_length,
        out_length,
        c,
        precision=DEFAULT_HIGH_PRECISION):
    """The coefficients of the DDT of a random function can be modeled like independent and identically distributed random variables following a specific Poisson distribution [JMC:DaeRij07]. This function returns this probability.

    This probability is identical for a random permutation and a random non-bijective function, and it corresponds to a Poisson distribution.

    Args:
        in_length (int): the input bit-length
        out_length (int): the output bit-length
        c (int): the DDT coefficient value whose probability we want
        precision (int): a facultative argument corresponding to the level of precision to use for the floating point arithmetic.

    Returns:
        RealNumber: The probability that a coefficient of the DDT of a function mapping `in_length` bits to `out_length` is equal to `coeff`.
    
    """
    m, n = in_length, out_length
    big_precision = RealField(precision)
    if c % 2 == 1:
        return 0
    else:
        d = c/2
        if m >= 5 and n-m <= m/4:
            k = 2**(m-n-1)
            return RealNumber(exp(-k) * k**d / factorial(d))
        else:
            return RealNumber(binomial(2**(m-1), d) * 2**(-n*d) * (1 - 2**-n)**(2**(m-1)-d))


def expected_differential_uniformity_distribution_permutation(
        in_length,
        out_length,
        vmin=2,
        vmax=28):
    """The coefficients of the DDT of a random function can be modeled like independent and identically distributed random variables following a specific Poisson distribution [JMC:DaeRij07]. As a consequence, we can estimate the value of the maximum coefficient.

    Args:
        in_length (int): the number of bits in the input
        out_length (int): the number of bits in the output
        vmin (int): the minimum value of the differential uniformity to consider
        vmax (int): the maximum value of the differential uniformity to consider

    Returns:
        dict: A dictionary `d` such that `d[i]` is the probability that an F_2 S-box mapping `in_length` bits to `out_length` has a differential uniformity equal to `i`.
    
    """
    return expected_maximum_distribution(
        in_length,
        out_length,
        ddt_coeff_probability,
        (2**in_length-1) * (2**out_length-1),
        vmin=vmin,
        vmax=vmax)

    
# !SUBSECTION! LAT


# !SUBSUBSECTION! Bijective case

def lat_coeff_probability_permutation(
        in_length,
        out_length,
        c,
        precision=DEFAULT_HIGH_PRECISION):
    """The non-trivial coefficients of the LAT of a permutation behave like samples from independent and identically distributed random variables following the same distribution, as explained in [JMC:DaeRij07] (where said distribution is also given). This distribution is not the same as in the case of non-bijective functions.

    If in_length != out_length, raises an error: a permutation cannot have a different input and output size.

    Args:
        in_length (int): the input bit-length
        out_length (int): the output bit-length
        c (int): the absolute value of the LAT coefficient value whose probability we want
        precision (int): a facultative argument corresponding to the level of precision to use for the floating point arithmetic.

    Returns:
        RealNumber: The probability that a coefficient of the LAT of a permutation operating on `in_length` bits has an absolute value equal to `c`.

    """
    m, n = in_length, out_length
    big_precision = RealField(precision)
    if m != n:
        raise "A permutation cannot map {} bits to {}!".format(m, n)
    if c % 4 != 0:
        return 0
    elif c == 0:
        return big_precision(Integer(binomial(2**(n-1), 2**(n-2)))**2) / Integer(binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return 2 * big_precision(Integer(binomial(2**(n-1), 2**(n-2) + c/2))**2) / Integer(binomial(2**n, 2**(n-1)))


def expected_linearity_distribution_permutation(
        in_length,
        out_length,
        vmin=8,
        vmax=100):
    """The coefficients of the LAT of a random permutation can be modeled like independent and identically distributed random variables following a specific Poisson distribution [JMC:DaeRij07]. As a consequence, we can estimate the value of the maximum absolute coefficient of the LAT (i.e., the linearity).

    Args:
        in_length (int): the number of bits in the input
        out_length (int): the number of bits in the output
        vmin (int): the minimum value of the linearity to consider
        vmax (int): the maximum value of the linearity to consider

    Returns:
        dict: A dictionary `d` such that `d[i]` is the probability that a permutation operating on `in_length` bits has a linearity uniformity equal to `i`.
    
    """
    return expected_maximum_distribution(
        in_length,
        out_length,
        lat_coeff_probability_permutation,
        (2**in_length-1) * (2**out_length-1),
        vmin=vmin,
        vmax=vmax)
    

# !SUBSUBSECTION! Non-bijective case

def lat_coeff_probability_function(
        in_length,
        out_length,
        c,
        precision=DEFAULT_HIGH_PRECISION):
    """The non-trivial coefficients of the LAT of a function behave like samples from independent and identically distributed random variables following the same distribution, as explained in [JMC:DaeRij07] (where said distribution is also given). This distribution is not the same as in the case of permutations.

    Args:
        in_length (int): the input bit-length
        out_length (int): the output bit-length
        c (int): the absolute value of the LAT coefficient value whose probability we want
        precision (int): a facultative argument corresponding to the level of precision to use for the floating point arithmetic.

    Returns:
        RealNumber: The probability that a coefficient of the LAT of a function mapping `in_length` bits to `out_length` bits has an absolute value equal to `c`.

    """
    m, n = in_length, out_length
    big_precision = RealField(precision)
    if m != n:
        raise "m (={}) should be equal to n (={})!".format(m, n)
    if c % 2 != 0:
        return 0
    if c == 0:
        return big_precision(2**(-2**n) * binomial(2**n, 2**(n-1)))
    else:
        c = c/2
        return 2 * big_precision(2**(-2**n) * binomial(2**n, 2**(n-1)+c))



def expected_linearity_distribution_function(
        in_length,
        out_length,
        vmin=8,
        vmax=100):
    """The coefficients of the LAT of a random function can be modeled like independent and identically distributed random variables following a specific Poisson distribution [JMC:DaeRij07]. As a consequence, we can estimate the value of the maximum absolute coefficient of the LAT (i.e., the linearity).

    Args:
        in_length (int): the number of bits in the input
        out_length (int): the number of bits in the output
        vmin (int): the minimum value of the linearity to consider
        vmax (int): the maximum value of the linearity to consider

    Returns:
        dict: A dictionary `d` such that `d[i]` is the probability that a permutation operating on `in_length` bits has a linearity uniformity equal to `i`.
    
    """
    return expected_maximum_distribution(
        in_length,
        out_length,
        lat_coeff_probability_permutation,
        (2**in_length-1) * 2**out_length,
        vmin=vmin,
        vmax=vmax)
    

# !SUBSECTION! BCT

# !TODO! Write proper docstring for the BCT anomalies functions

def bct_coeff_probability(
        in_length,
        out_length,
        c,
        precision=DEFAULT_HIGH_PRECISION):
    """The non-trivial coefficients of the BCT of a permutation behave like samples from independent and identically distributed random variables following the same distribution, as explained in [AC:BonPerTia19] (where said distribution is also given). 

    If in_length != out_length, raises an error: the BCT is not defined in this case.

    Args:
        in_length (int): the input bit-length
        out_length (int): the output bit-length
        c (int): the absolute value of the BCT coefficient value whose probability we want
        precision (int): a facultative argument corresponding to the level of precision to use for the floating point arithmetic.

    Returns:
        RealNumber: The probability that a coefficient of the BCT of a permutation operating on `in_length` bits has an absolute value equal to `c`.

    """
    m, n = in_length, out_length
    big_precision = RealField(precision)
    if m != n:
        raise "the BCT is only defined when in_length == out_length"
    if c % 2 == 1:
        return RealNumber(0.0)
    B = big_precision(2**(n-1))
    A = big_precision(2**(2*n-2)-2**(n-1))
    p = big_precision(1/(2**n-1))
    q = p**2
    d = int(c/2)
    result = big_precision(0.0)
    base = big_precision(2**n-1)**(B + 2*A)
    for j1, j2 in itertools.product(range(0, d+1), range(0, d+1)):
        if 2*j1 + 4*j2 == c:
            added = big_precision(binomial(B, j1)) * big_precision((2**n-2)**(B-j1)) * big_precision(binomial(A, j2)) * big_precision((2**(2*n)-2**(n+1))**(A-j2))
            result += added / base
    return result


def expected_boomerang_uniformity_distribution_permutation(
        in_length,
        out_length,
        vmin=2,
        vmax=28):
    """The coefficients of the BCT of a random function can be modeled like independent and identically distributed random variables following a specific Poisson distribution [AC:BonPerTia19]. As a consequence, we can estimate the value of the maximum coefficient.

    Args:
        in_length (int): the number of bits in the input
        out_length (int): the number of bits in the output
        vmin (int): the minimum value of the differential uniformity to consider
        vmax (int): the maximum value of the differential uniformity to consider

    Returns:
        dict: A dictionary `d` such that `d[i]` is the probability that an F_2 S-box mapping `in_length` bits to `out_length` has a boomerang uniformity equal to `i`.
    
    """
    return expected_maximum_distribution(
        in_length,
        out_length,
        bct_coeff_probability,
        (2**in_length-1) * (2**out_length-1),
        vmin=vmin,
        vmax=vmax)
    
# !SECTION! Aggregated information from the tables
# !TODO! Write proper docstring for the high level anomalies functions


# !SUBSECTION! Helper functions 

def probability_of_max_and_occurrences(
        in_length,
        out_length,
        v_max,
        occurrences,
        proba_func,
        precision=DEFAULT_HIGH_PRECISION):
    """Returns the logarithm in base 2 of the probability that
    $(2^m-1)(2^n-1)$ trials of an experiment yielding output c with
    probability proba_func(m, n, c) will have a result equal to
    `v_max` at most `occurrences` times and be strictly smaller on all
    other trials.

    """
    try:
        m, n = in_length, out_length
        big_precision = RealField(precision)
        p_strictly_smaller = sum(big_precision(proba_func(m, n, i, precision=precision)) for i in range(0, v_max))
        p_equal = big_precision(proba_func(m, n, v_max, precision=precision))
        result = big_precision(0)
        n_trials = (2**m-1) * (2**n-1)
        for occ in reversed(range(0, occurrences+1)):
            added = big_precision(binomial(n_trials, occ)) * (p_strictly_smaller)**((n_trials - occ)) * p_equal**(occ)
            result += added
        if imag_part(result) != 0: # can happen sometimes, I (Léo) have no idea why
            raise Exception("somehow, an imaginary number returned")
        return result
    except:
        print("failure in sum: occ={}, result={}, added={}".format(occ, result, added))
        return 0



def get_proba_func(s, table):
    if table == "DDT":
        return ddt_coeff_probability
    elif table == "LAT":
        if is_permutation(s):
            return lat_coeff_probability_permutation
        else:
            return lat_coeff_probability_function
    elif table == "BCT":
        return bct_coeff_probability
    else:
        raise Exception("unknown table name: " + table_name)
    

# !SUBSECTION! Anomaly computations

def table_anomaly(
        s,
        table,
        spec=None,
        precision=DEFAULT_HIGH_PRECISION):
    """Computes the positive anomaly (in the sense of [AC:BonPerTia19]) of the S-boxable object `s` that corresponds to its DDT, LAT or BCT.

    Args:
        s: the S_boxable object you want the anomaly of.
        table: the name of the table for which the anomaly must be computed. Must be either "DDT", "LAT" or "BCT".

    Returns:
        RealNumber: the positive anomaly associated to the given table for the S-box corresponding to `s`.
    
    """
    if table not in ["DDT", "LAT", "BCT"]:
        raise "The table-based anomaly is defined for the LAT, DDT and BCT. table={} is unknown.".format(table)
    else:
        sb = Sb(s)
        proba_func = get_proba_func(sb, table)
        if table == "DDT":
            if spec == None:
                spec = differential_spectrum(s)
        elif table == "LAT":
            if spec == None:
                spec = walsh_spectrum(s)
        else:
            if spec == None:
                spec = boomerang_spectrum(s)
        v_max = 0
        occurrences = 0
        for k in spec.keys():
            if abs(k) == v_max:
                occurrences += spec[k]
            elif abs(k) > v_max:
                v_max = abs(k)
                occurrences = spec[k]
        big_precision = RealField(precision)
        p = big_precision(
            probability_of_max_and_occurrences(
                sb.get_input_length(),
                sb.get_output_length(),
                v_max,
                occurrences,
                proba_func,
                precision=precision)
        )
        return -p.log2()

    
def table_negative_anomaly(
        s,
        table,
        spec=None,
        precision=DEFAULT_HIGH_PRECISION):
    """Computes the negative anomaly (in the sense of [AC:BonPerTia19]) of the S_box `s` that corresponds to its DDT, LAT or BCT.

    Args:
        s: the S_boxable object you want the anomaly of.
        table: the name of the table for which the anomaly must be computed. Must be either "DDT", "LAT" or "BCT".

    Returns:
        RealNumber: the negative anomaly associated to the given table for the S-box corresponding to `s`.
    
    """
    if table not in ["DDT", "LAT", "BCT"]:
        raise "The table-based negative anomaly is defined for the LAT, DDT and BCT. table={} is unknown.".format(table)
    else:
        sb = Sb(s)
        proba_func = get_proba_func(sb, table)
        m, n = sb.get_input_length(), sb.get_output_length()
        if table == "DDT":
            if spec == None:
                spec = differential_spectrum(s)
        elif table == "LAT":
            if spec == None:
                spec = walsh_spectrum(s)
        else:
            if spec == None:
                spec = boomerang_spectrum(s)
        v_max = 0
        occurrences = 0
        for k in spec.keys():
            if abs(k) == v_max:
                occurrences += spec[k]
            elif abs(k) > v_max:
                v_max = abs(k)
                occurrences = spec[k]
        p_anomaly = probability_of_max_and_occurrences(
            m,
            n,
            v_max,
            occurrences,
            proba_func,
            precision=precision
        )
        p_equal = proba_func(m, n, v_max, precision=precision)
        p_strictly_smaller = sum(proba_func(m, n, i, precision=precision) for i in range(0, v_max))
        card = (2**n-1)**2
        big_precision = RealField(precision)
        p_precise_equal = big_precision(binomial(card, occurrences))*p_equal**occurrences*p_strictly_smaller**(card-occurrences)
        # p_precise_equal = min(p_precise_equal, big_precision(1.0))
        # p_anomaly       = min(p_anomaly, big_precision(1.0))
        return -big_precision(p_precise_equal-p_anomaly+1).log2()

