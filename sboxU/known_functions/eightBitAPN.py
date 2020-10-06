"""Contains many 8-bit APN functions"""

from sage.all import GF, PolynomialRing

# global variables of the module
N = 8
F = GF(2**N, name="a")
g = F.gen()
POLY_RING = PolynomialRing(F, "X")
X = POLY_RING.gen()


def poly_to_lut(p):
    s = []
    for x_i in xrange(0, 2**N):
        y = (p(F.fetch_int(x_i))).integer_representation()
        s.append(y)
    return s


def all_quadratic_polynomials():
    """All the functions in Table 9 of
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.215.5432&rep=rep1&type=pdf

    """
    return [
        poly_to_lut(X**3),
        poly_to_lut(X**9),
        poly_to_lut(X**3 + sum(X**(9*2**i) for i in xrange(0, N))),
        poly_to_lut(X**9 + sum(X**(3*2**i) for i in xrange(0, N))),
        poly_to_lut(X**3 + g**245*X**33 + g**183*X**66 + g**21*X**144),
        poly_to_lut(X**3 + g**65*X**18 + g**120*X**66 + g**135*X**144),

        poly_to_lut(g**188*X**192 + g**129*X**144 + g**172*X**132 + g**138*X**129 + g**74*X**96 + g**244*X**72 + g**22*X**66 + g**178*X**48 + g**150*X**36 + g**146*X**33 + g**6*X**24 + g**60*X**18 + g**80*X**12 + g**140*X**9 + g**221*X**6 + g**19*X**3),
        poly_to_lut(g**37*X**192 + g**110*X**144 + g**40*X**132 + g**53*X**129 + g**239*X**96 + g**235*X**72 + g**126*X**66 + g**215*X**48 + g**96*X**36 + g**29*X**33 + g**19*X**24 + g**14*X**18 + g**139*X**12 + g**230*X**9 + g**234*X**6 + g**228*X**3),
        poly_to_lut(g**242*X**192 + g**100*X**144 + g**66*X**132 + g**230*X**129 + g**202*X**96 + g**156*X**72 + g**254*X**66 + g**18*X**48 + g**44*X**36 + g**95*X**33 + g**100*X**24 + g**245*X**18 + g**174*X**12 + g**175*X**9 + g**247*X**6 + g**166*X**3),
        poly_to_lut(g**100*X**192 + g**83*X**144 + g**153*X**132 + g**65*X**129 + g**174*X**96 + g**136*X**72 + g**46*X**66 + g**55*X**48 + g**224*X**36 + g**180*X**33 + g**179*X**24 + g**226*X**18 + g**54*X**12 + g**168*X**9 + g**89*X**6 + g**56*X**3),
        poly_to_lut(g**77*X**192 + g**133*X**144 + g**47*X**132 + g**229*X**129 + g**23*X**96 + g**242*X**72 + g**242*X**66 + g**245*X**48 + g**212*X**36 + g**231*X**33 + g**174*X**24 + g**216*X**18 + g**96*X**12 + g**253*X**9 + g**154*X**6 + g**71*X**3),
        poly_to_lut(g**220*X**192 + g**94*X**144 + g**70*X**132 + g**159*X**129 + g**145*X**96 + g**160*X**72 + g**74*X**66 + g**184*X**48 + g**119*X**36 + g**106*X**33 + g**253*X**24 + g*X**18 + g**90*X**12 + g**169*X**9 + g**118*X**6 +   + g**187*X**3),
        poly_to_lut(g**98*X**192 + g**225*X**144 + g**111*X**132 + g**238*X**129 + g**182*X**96 + g**125*X**72 + g**196*X**66 + g**219*X**48 + g**189*X**36 + g**199*X**33 + g**181*X**24 + g**110*X**18 + g**19*X**12 + g**175*X**9 + g**133*X**6 + g**47*X**3),
        poly_to_lut(g**236*X**192 + g**212*X**160 + g**153*X**144 + g**185*X**136 + g**3*X**132 + g**89*X**130 + g**189*X**129 + g**182*X**96 + g**105*X**80 + g**232*X**72 + g**219*X**68 + g**145*X**66 + g**171*X**65 + g**107*X**48 + g**179*X**40 + g**227*X**36 + g**236*X**34 + g**189*X**33 + g**162*X**24 + g**216*X**20 + g**162*X**18 + g**117*X**17 + g**56*X**12 + g**107*X**10 + g**236*X**9 + g**253*X**6 + g**180*X**5 + g**18*X**3),
        poly_to_lut(g**27*X**192 + g**167*X**144 + g**26*X**132 + g**231*X**129 + g**139*X**96 + g**30*X**72 + g**139*X**66 + g**203*X**48 + g**36*X**36 + g**210*X**33 + g**195*X**24 + g**12*X**18 + g**43*X**12 + g**97*X**9 + g**61*X**6 + g**39*X**3),
        poly_to_lut(g**6*X**192 + g**85*X**144 + g**251*X**132 + g**215*X**129 + g**229*X**96 + g**195*X**72 + g**152*X**66 + g**173*X**48 + g**209*X**36 + g**165*X**33 + g**213*X**24 + g**214*X**18 + g**158*X**12 + g**146*X**9 + X**6 + g**50*X**3),
        poly_to_lut(g**164*X**192 + g**224*X**144 + g**59*X**132 + g**124*X**129 + g**207*X**96 + g**211*X**72 + g**5*X**66 + g**26*X**48 + g**20*X**36 + g**101*X**33 + g**175*X**24 + g**241*X**18 + X**12 + g**15*X**9 + g**217*X**6 + g**212*X**3),
        
        poly_to_lut(X**3 + X**17 + g**16*(X**18 + X**33) + g**15*X**48),
        poly_to_lut(X**3 + g**24*X**6 + g**182*X**132 + g**67*X**192),
        poly_to_lut(X**3 + X**6 + X**68 + X**80 + X**132 + X**160),
        poly_to_lut(X**3 + X**5 + X**18 + X**40 + X**66),
        poly_to_lut(X**3 + X**12 + X**40 + X**66 + X**130),
    ]


def all_QAMs():
    """All the functions found using the QAM method, see
    https://link.springer.com/article/10.1007/s10623-014-9955-3
    
    """
    from allQam import all_funcs
    return all_funcs


def all_BeiLea():
    """All the functions found by Beierle and Leander in 2020, availabe
    online at: https://zenodo.org/record/4030734

    """
    from BeierleLeander import all_funcs
    return all_funcs


def all_quadratics():
    """Returns all known quadratic APN functions operating on 8 bits. They
    are all in distinct CCZ-classes.

    """
    return all_QAMs() + all_BeiLea()


def all_non_quadratics():
    return [
        poly_to_lut(X**57),
    ]
