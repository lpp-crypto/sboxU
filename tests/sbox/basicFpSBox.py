#!/usr/bin/env python
import sys
from sage.all import *
from sboxU import *
from sage.all import GF, PolynomialRing

from sboxU.core.sbox.cython_functions import S_box_fp

p = 3

Fp = GF(p)

raw_lut = [[0, 1], [1, 0], [0, 2], [0, 0], [2, 0], [2, 2], [1, 2], [2, 1], [1, 1]]

lut = [[Fp(x) for x in row] for row in raw_lut]


def fp_eq(a, b):
    """Return True if two FpWord lists are equal entry by entry."""

    return list(a) == list(b)


def fp_lut_eq(lut_a, lut_b):
    """Return True if two full LUTs are equal."""

    return all(fp_eq(a, b) for a, b in zip(lut_a, lut_b))


def main_test():
    with Experiment("S-boxes over F_p"):
        section("Construction")
        subsection("From a lookup table")
        # --- {
        u = get_sbox(lut)
        if isinstance(u, S_box_fp):
            success("get_sbox(lut) correctly returned an S_box_fp instance")
        else:
            fail("expected S_box_fp, got {}".format(type(u)))
        # --- }
        subsection("From an S_box_fp")
        sprime = get_sbox(u, name="rename of s")
        sprime_name = sprime.get_name()
        if sprime_name == b"rename of s":
            success(
                "Passing a name argument in get_sbox when input is of type S_box_fp creates an sbox with a new name"
            )
        else:
            fail(
                f"something went wrong, the name of the new sbox is not as expected : name = {sprime_name}"
            )

        subsection("From multivariate polynomials")
        # --- {
        R = PolynomialRing(Fp, 2, "x")
        x1, x2 = R.gens()
        v = get_sbox([x1 * x2, x1 + x2])
        if isinstance(v, S_box_fp):
            success("get_sbox(polynomials) correctly returned an S_box_fp instance")
        else:
            fail("expected S_box_fp from polynomial input, got {}".format(type(v)))
        # --- }
        # --- {
        if fp_eq(v[[0, 0]], [0, 0]):
            success("v[(0,0)] = (0,0) as expected")
        else:
            fail("v[(0,0)] = {} instead of (0,0)".format(list(v[[0, 0]])))
        # --- }
        # --- {
        if fp_eq(v[[1, 2]], [2, 0]):
            success("v[(1,2)] = (2,0) as expected")
        else:
            fail("v[(1,2)] = {} instead of (2,0)".format(list(v[[1, 2]])))
        # --- }
        subsection("From a univariate polynomial")
        # --- {
        P3 = PolynomialRing(Fp, "x")
        xp = P3.gen()
        sq_uni = get_sbox(xp**2)
        if isinstance(sq_uni, S_box_fp):
            success("get_sbox(x**2 over GF(3)) returned an S_box_fp instance")
        else:
            fail(
                "expected S_box_fp from univariate polynomial, got {}".format(
                    type(sq_uni)
                )
            )
        # --- }
        # --- {
        expected_sq = [[0], [1], [1]]
        got_sq = [list(y) for y in sq_uni.get_lut()]
        if got_sq == expected_sq:
            success("squaring map x^2 over GF(3) has LUT [[0],[1],[1]]")
        else:
            fail("squaring map LUT is {}, expected {}".format(got_sq, expected_sq))
        # --- }
        # --- {
        id_uni = get_sbox(xp)
        id_sp = [list(v) for v in id_uni.get_input_space()]
        id_lut = [list(y) for y in id_uni.get_lut()]
        if id_lut == id_sp:
            success("identity polynomial x over GF(3) gives the identity S-box")
        else:
            fail("identity polynomial LUT is {}, expected {}".format(id_lut, id_sp))
        # --- }
        # --- {
        K9 = GF(9, "a")
        P9 = PolynomialRing(K9, "z")
        z = P9.gen()
        id9_uni = get_sbox(z)
        sp9 = [list(v) for v in id9_uni.get_input_space()]
        lut9 = [list(y) for y in id9_uni.get_lut()]
        if lut9 == sp9 and id9_uni.get_input_size() == 2:
            success("identity polynomial z over GF(9) gives identity S-box on F_3^2")
        else:
            fail(
                "GF(9) identity polynomial: input_size={}, lut==input_space: {}".format(
                    id9_uni.get_input_size(), lut9 == sp9
                )
            )
        # --- }
        subsection("Bytes round-trip")
        # --- {
        b = u.to_bytes()
        if b[0] == 0:
            success("leading byte of to_bytes() is 0x00 as required by the Fp format")
        else:
            fail("leading byte of to_bytes() is {}, expected 0".format(b[0]))
        # --- }
        # --- {
        u_reconstructed = get_sbox(b)
        if isinstance(u_reconstructed, S_box_fp):
            success("get_sbox(to_bytes()) returns an S_box_fp")
        else:
            fail(
                "get_sbox(to_bytes()) returned {}, expected S_box_fp".format(
                    type(u_reconstructed)
                )
            )
        # --- }
        # --- {
        if u_reconstructed == u:
            success("round-trip to_bytes -> get_sbox is the identity")
        else:
            fail("round-trip failed: reconstructed S-box differs from original")
        # --- }
        section("Basic properties")
        subsection("Characteristic and sizes")
        # --- {
        if u.get_p() == 3:
            success("get_p() returns 3")
        else:
            fail("get_p() returned {}, expected 3".format(u.get_p()))
        # --- }
        # --- {
        if u.get_input_size() == 2:
            success("get_input_size() returns 2")
        else:
            fail("get_input_size() returned {}, expected 2".format(u.get_input_size()))
        # --- }
        # --- {
        if u.get_output_size() == 2:
            success("get_output_size() returns 2")
        else:
            fail(
                "get_output_size() returned {}, expected 2".format(u.get_output_size())
            )
        # --- }
        # --- {
        if u.get_input_space_size() == 9:
            success("get_input_space_size() returns p^t = 9")
        else:
            fail(
                "get_input_space_size() returned {}, expected 9".format(
                    u.get_input_space_size()
                )
            )
        # --- }
        # --- {
        if u.get_output_space_size() == 9:
            success("get_output_space_size() returns p^u = 9")
        else:
            fail(
                "get_output_space_size() returned {}, expected 9".format(
                    u.get_output_space_size()
                )
            )
        # --- }
        # --- {
        if len(u) == 9:
            success("len(u) returns 9")
        else:
            fail("len(u) returned {}, expected 9".format(len(u)))
        # --- }
        subsection("Input and output spaces")
        # --- {
        expected_input_space = [
            [0, 0],
            [1, 0],
            [2, 0],
            [0, 1],
            [1, 1],
            [2, 1],
            [0, 2],
            [1, 2],
            [2, 2],
        ]
        actual = [list(x) for x in u.get_input_space()]
        if actual == expected_input_space:
            success(
                "get_input_space() enumerates all elements in correct (little-endian) order"
            )
        else:
            fail(
                "get_input_space() returned {}, expected {}".format(
                    actual, expected_input_space
                )
            )
        # --- }
        # --- {
        if [list(x) for x in u.get_output_space()] == expected_input_space:
            success("get_output_space() matches the input space (same dimensions)")
        else:
            fail("get_output_space() has unexpected content")
        # --- }
        subsection("Lookup table access")
        # --- {
        lut_u = [list(x) for x in u.get_lut()]
        if lut_u == raw_lut:
            success("get_lut() returns the original lookup table")
        else:
            fail("get_lut() returned {}, expected {}".format(lut_u, raw_lut))
        # --- }
        subsection("String representation")
        # --- {
        expected_str = "[[0,1],[1,0],[0,2],[0,0],[2,0],[2,2],[1,2],[2,1],[1,1]]"
        if str(u) == expected_str:
            success("str(u) matches the expected representation")
        else:
            fail("str(u) = '{}', expected '{}'".format(str(u), expected_str))
        # --- }
        section("Evaluation")
        subsection("Direct lookup with __getitem__")
        # --- {
        checks = [
            ([0, 0], [0, 1]),
            ([1, 0], [1, 0]),
            ([2, 0], [0, 2]),
            ([0, 1], [0, 0]),
            ([1, 1], [2, 0]),
            ([2, 1], [2, 2]),
            ([0, 2], [1, 2]),
            ([1, 2], [2, 1]),
            ([2, 2], [1, 1]),
        ]
        all_ok = True
        for x, expected in checks:
            if not fp_eq(u[x], expected):
                fail("u[{}] = {}, expected {}".format(x, list(u[x]), expected))
                all_ok = False
        if all_ok:
            success("u[x] is correct for all 9 inputs")
        # --- }
        subsection("Iteration")
        # --- {
        iterated = [list(y) for y in u]
        if iterated == raw_lut:
            success("iterating over u yields all outputs in input-space order")
        else:
            fail("iteration gave {}, expected {}".format(iterated, raw_lut))
        # --- }
        section("Equality, inequality and hash")
        # --- {
        u2 = get_sbox(lut)
        if u == u2:
            success("two S-boxes built from the same LUT are equal")
        else:
            fail("u == u2 returned False, expected True")
        # --- }
        # --- {
        if u != v:
            success("u != v when the two S-boxes differ")
        else:
            fail("u != v returned False for distinct S-boxes")
        # --- }
        # --- {
        if hash(u) == hash(u2):
            success("equal S-boxes have the same hash")
        else:
            fail("equal S-boxes have different hashes")
        # --- }
        section("Arithmetic operators")
        subsection("Addition")
        # --- {
        uu = u + u
        expected_double = [[(2 * a) % 3 for a in row] for row in raw_lut]
        got_double = [list(y) for y in uu]
        if got_double == expected_double:
            success("u + u doubles every output coordinate mod 3")
        else:
            fail("u + u gave {}, expected {}".format(got_double, expected_double))
        # --- }
        subsection("Subtraction")
        # --- {
        u_minus_u = u - u
        zero_lut = [[0, 0]] * 9
        got_zero = [list(y) for y in u_minus_u]
        if got_zero == zero_lut:
            success("u - u is the zero function")
        else:
            fail("u - u gave {}, expected all zeros".format(got_zero))
        # --- }
        # --- {
        upv_minus_v = (u + v) - v
        if upv_minus_v == u:
            success("(u + v) - v == u")
        else:
            fail("(u + v) - v != u: got {}".format([list(y) for y in upv_minus_v]))
        # --- }
        subsection("Composition")
        # --- {
        uu_comp = u * u
        expected_comp = [list(u[list(row)]) for row in raw_lut]
        got_comp = [list(y) for y in uu_comp]
        if got_comp == expected_comp:
            success("(u * u)[x] == u(u(x)) for all x")
        else:
            fail("composition gave {}, expected {}".format(got_comp, expected_comp))
        # --- }
        subsection("Power operator")
        # --- {
        id_u = u**0
        identity_lut = [list(x) for x in u.get_input_space()]
        got_id = [list(y) for y in id_u]
        if got_id == identity_lut:
            success("u ** 0 is the identity")
        else:
            fail("u ** 0 gave {}, expected identity {}".format(got_id, identity_lut))
        # --- }
        # --- {
        if (u**2) == (u * u):
            success("u ** 2 == u * u")
        else:
            fail("u ** 2 != u * u")
        # --- }
        # --- {
        if (u**-1) == u.inverse():
            success("u ** -1 == u.inverse()")
        else:
            fail("u ** -1 != u.inverse()")
        # --- }
        section("Inversion")
        subsection("Identity S-box")
        # --- {
        id_sb = S_box_fp.identity_S_box(2, 3)
        if fp_lut_eq(id_sb.get_lut(), u.get_input_space()):
            success("identity_S_box(2,3) maps every input to itself")
        else:
            fail(
                "identity_S_box(2,3) has wrong LUT: {}".format(
                    [list(x) for x in id_sb.get_lut()]
                )
            )
        # --- }
        subsection("Invertibility test")
        # --- {
        if u.is_invertible():
            success("u.is_invertible() returns True for the bijection u")
        else:
            fail("u.is_invertible() returned False for a known bijection")
        # --- }
        subsection("Inverse function")
        # --- {
        u_inv = u.inverse()
        id_sb = S_box_fp.identity_S_box(2, 3)
        # --- }
        # --- {
        if (u * u_inv) == id_sb:
            success("u * u.inverse() is the identity")
        else:
            fail("u * u.inverse() is not the identity")
        # --- }
        # --- {
        if (u_inv * u) == id_sb:
            success("u.inverse() * u is the identity")
        else:
            fail("u.inverse() * u is not the identity")
        # --- }
        # --- {
        if fp_eq(u_inv[[0, 0]], [0, 1]):
            success("u_inv[(0,0)] = (0,1) as expected from the LUT")
        else:
            fail("u_inv[(0,0)] = {}, expected (0,1)".format(list(u_inv[[0, 0]])))
        # --- }
        section("Coordinate functions and components")
        subsection("Coordinate functions")
        # --- {
        coord_0 = u.coordinate(0)
        expected_coord_0 = [[row[0]] for row in raw_lut]
        got_coord_0 = [list(y) for y in coord_0]
        if got_coord_0 == expected_coord_0:
            success("coordinate(0) correctly extracts the first output coordinate")
        else:
            fail(
                "coordinate(0) gave {}, expected {}".format(
                    got_coord_0, expected_coord_0
                )
            )
        # --- }
        # --- {
        coord_1 = u.coordinate(1)
        expected_coord_1 = [[row[1]] for row in raw_lut]
        got_coord_1 = [list(y) for y in coord_1]
        if got_coord_1 == expected_coord_1:
            success("coordinate(1) correctly extracts the second output coordinate")
        else:
            fail(
                "coordinate(1) gave {}, expected {}".format(
                    got_coord_1, expected_coord_1
                )
            )
        # --- }
        subsection("Component functions")
        # --- {
        comp_10 = u.component([1, 0])
        if comp_10 == coord_0:
            success("component([1,0]) equals coordinate(0)")
        else:
            fail(
                "component([1,0]) differs from coordinate(0): {} vs {}".format(
                    [list(y) for y in comp_10], [list(y) for y in coord_0]
                )
            )
        # --- }
        # --- {
        comp_01 = u.component([0, 1])
        if comp_01 == coord_1:
            success("component([0,1]) equals coordinate(1)")
        else:
            fail("component([0,1]) differs from coordinate(1)")
        # --- }
        # --- {
        expected_comp_11 = [[(row[0] + row[1]) % 3] for row in raw_lut]
        got_comp_11 = [list(y) for y in u.component([1, 1])]
        if got_comp_11 == expected_comp_11:
            success("component([1,1]) computes the sum of coordinates mod 3")
        else:
            fail(
                "component([1,1]) gave {}, expected {}".format(
                    got_comp_11, expected_comp_11
                )
            )
        # --- }
        # --- {
        expected_comp_21 = [[(2 * row[0] + row[1]) % 3] for row in raw_lut]
        got_comp_21 = [list(y) for y in u.component([2, 1])]
        if got_comp_21 == expected_comp_21:
            success("component([2,1]) computes 2*s0 + s1 mod 3")
        else:
            fail(
                "component([2,1]) gave {}, expected {}".format(
                    got_comp_21, expected_comp_21
                )
            )
        # --- }
        # --- {
        comp_00 = u.component([0, 0])
        if [list(y) for y in comp_00] == [[0]] * 9:
            success("component([0,0]) is the zero function")
        else:
            fail(
                "component([0,0]) is not the zero function: {}".format(
                    [list(y) for y in comp_00]
                )
            )
        # --- }
        section("Derivative")
        # --- {
        zero_delta = [0, 0]
        d_zero = u.derivative(zero_delta)
        if [list(y) for y in d_zero] == [[0, 0]] * 9:
            success("derivative in direction 0 is the zero function")
        else:
            fail(
                "derivative in direction 0 gave {}, expected all zeros".format(
                    [list(y) for y in d_zero]
                )
            )
        # --- }
        # --- {
        delta = [1, 0]
        input_sp = [list(x) for x in u.get_input_space()]
        input_size = u.get_input_size()
        n = u.get_input_space_size()
        powers = u.get_powers_in()
        expected_deriv = []
        for i in range(n):
            x = input_sp[i]
            x_plus_delta = [Fp(x[j]) + Fp(delta[j]) for j in range(input_size)]
            idx_xd = S_box_fp.vec_to_int(x_plus_delta, powers)
            out_xd = list(u.get_lut()[idx_xd])
            out_x = list(u.get_lut()[i])
            expected_deriv.append(
                [Fp(out_xd[j]) - Fp(out_x[j]) for j in range(u.get_output_size())]
            )
        got_deriv = [list(y) for y in u.derivative(delta)]
        if got_deriv == expected_deriv:
            success("derivative in direction (1,0) matches the expected values")
        else:
            fail(
                "derivative in direction (1,0) gave {}, expected {}".format(
                    got_deriv, expected_deriv
                )
            )
        # --- }
        section("Serialization")
        # --- {
        import struct

        bs = u.to_bytes()
        # byte 0: format marker
        # bytes 1..8: p (8-byte LE)
        # bytes 9..16: output_size (8-byte LE)
        # bytes 17..24: n = p^t (8-byte LE)
        # byte 25: bytes_per_value
        p_encoded = struct.unpack_from("<q", bs, 1)[0]
        out_encoded = struct.unpack_from("<q", bs, 9)[0]
        n_encoded = struct.unpack_from("<q", bs, 17)[0]
        bpv = bs[25]
        # --- }
        # --- {
        if p_encoded == 3:
            success("p is correctly encoded as 3 in the header")
        else:
            fail("p encoded as {}, expected 3".format(p_encoded))
        # --- }
        # --- {
        if out_encoded == 2:
            success("output_size is correctly encoded as 2 in the header")
        else:
            fail("output_size encoded as {}, expected 2".format(out_encoded))
        # --- }
        # --- {
        if n_encoded == 9:
            success("n = p^t = 9 is correctly encoded in the header")
        else:
            fail("n encoded as {}, expected 9".format(n_encoded))
        # --- }
        # --- {
        if bpv == 1:
            success("bytes_per_value is 1 (p=3 fits in a single byte)")
        else:
            fail("bytes_per_value is {}, expected 1".format(bpv))
        # --- }
        # --- {
        expected_len = 1 + 3 * 8 + 1 + 9 * 2 * 1
        if len(bs) == expected_len:
            success("byte length of to_bytes() is {} as expected".format(expected_len))
        else:
            fail("byte length is {}, expected {}".format(len(bs), expected_len))
        # --- }
        section("S-boxes over other primes")
        subsection("p = 5, bijection")
        # --- {
        p5 = 5
        Fp5 = GF(p5)
        raw_lut_f = [[0], [2], [4], [1], [3]]
        lut_f = [[Fp5(x) for x in row] for row in raw_lut_f]
        f = get_sbox(lut_f)
        if isinstance(f, S_box_fp):
            success("p=5 S-box is correctly an S_box_fp instance")
        else:
            fail("p=5 S-box: expected S_box_fp, got {}".format(type(f)))
        # --- }
        # --- {
        if (
            f.get_p() == 5
            and f.get_input_size() == 1
            and f.get_output_size() == 1
            and len(f) == 5
        ):
            success("p=5 S-box has correct characteristic, dimensions and length")
        else:
            fail(
                "p=5 S-box: wrong properties p={}, t={}, u={}, len={}".format(
                    f.get_p(), f.get_input_size(), f.get_output_size(), len(f)
                )
            )
        # --- }
        # --- {
        if (
            fp_eq(f[[0]], [0])
            and fp_eq(f[[1]], [2])
            and fp_eq(f[[3]], [1])
            and fp_eq(f[[4]], [3])
        ):
            success("f[x] = 2x mod 5 is evaluated correctly")
        else:
            fail("f evaluation incorrect")
        # --- }
        # --- {
        if f.is_invertible():
            success("multiplication by 2 in F_5 is correctly identified as a bijection")
        else:
            fail("multiplication by 2 in F_5 should be a bijection")
        # --- }
        # --- {
        ff = f + f
        expected_4x = [[(4 * x) % 5] for x in range(5)]
        if [list(y) for y in ff] == expected_4x:
            success("f + f = 4x mod 5 for all x in F_5")
        else:
            fail(
                "f + f gave {}, expected {}".format([list(y) for y in ff], expected_4x)
            )
        # --- }
        # --- {
        if [list(y) for y in (f - f)] == [[0]] * 5:
            success("f - f is the zero function over F_5")
        else:
            fail("f - f is not zero over F_5")
        # --- }
        # --- {
        f_rt = get_sbox(f.to_bytes())
        if isinstance(f_rt, S_box_fp) and f_rt == f:
            success("p=5 S-box survives to_bytes -> get_sbox round-trip")
        else:
            fail("p=5 round-trip failed")
        # --- }
        subsection("p = 5, non-bijection")
        # --- {
        raw_lut_g = [[0], [1], [4], [4], [1]]
        lut_g = [[Fp5(x) for x in row] for row in raw_lut_g]
        g = get_sbox(lut_g)
        if not g.is_invertible():
            success("x^2 mod 5 is correctly identified as non-invertible")
        else:
            fail("x^2 mod 5 should not be invertible")
        # --- }
        subsection("p = 7, round-trip")
        # --- {
        p7 = 7
        Fp7 = GF(p7)
        raw_lut_h7 = [[(x + 1) % p7] for x in range(p7)]
        lut_h7 = [[Fp7(x) for x in row] for row in raw_lut_h7]
        h7 = get_sbox(lut_h7)
        h7_rt = get_sbox(h7.to_bytes())
        if isinstance(h7_rt, S_box_fp) and h7_rt == h7:
            success("p=7 S-box survives to_bytes -> get_sbox round-trip")
        else:
            fail("p=7 round-trip failed")
        # --- }
        # --- {
        bs7 = h7.to_bytes()
        p7_encoded = struct.unpack_from("<q", bs7, 1)[0]
        if p7_encoded == 7:
            success("p=7 is correctly encoded in the bytes header")
        else:
            fail("p=7 header encodes {}, expected 7".format(p7_encoded))
        # --- }
        section("Non-square S-boxes")
        subsection("Injective: F_3^1 to F_3^2")
        # --- {
        raw_lut_phi = [[0, 0], [1, 2], [2, 1]]
        lut_phi = [[Fp(x) for x in row] for row in raw_lut_phi]
        phi = get_sbox(lut_phi)
        if phi.get_input_size() == 1 and phi.get_output_size() == 2:
            success("phi: F_3^1 -> F_3^2 has correct dimensions (t=1, u=2)")
        else:
            fail(
                "phi dimensions wrong: t={}, u={}".format(
                    phi.get_input_size(), phi.get_output_size()
                )
            )
        # --- }
        # --- {
        if (
            len(phi) == 3
            and phi.get_input_space_size() == 3
            and phi.get_output_space_size() == 9
        ):
            success("phi has 3 inputs and 9 possible outputs")
        else:
            fail(
                "phi wrong space sizes: len={}, in={}, out={}".format(
                    len(phi), phi.get_input_space_size(), phi.get_output_space_size()
                )
            )
        # --- }
        # --- {
        if (
            fp_eq(phi[[0]], [0, 0])
            and fp_eq(phi[[1]], [1, 2])
            and fp_eq(phi[[2]], [2, 1])
        ):
            success("phi[x] = (x, 2x mod 3) evaluated correctly on all 3 inputs")
        else:
            fail(
                "phi evaluation failed: phi(0)={}, phi(1)={}, phi(2)={}".format(
                    list(phi[[0]]), list(phi[[1]]), list(phi[[2]])
                )
            )
        # --- }
        # --- {
        if not phi.is_invertible():
            success("phi is correctly non-invertible (input_size != output_size)")
        else:
            fail("phi should not be invertible since input_size != output_size")
        # --- }
        # --- {
        phi_rt = get_sbox(phi.to_bytes())
        if isinstance(phi_rt, S_box_fp) and phi_rt == phi:
            success("non-square S-box (t=1, u=2) survives bytes round-trip")
        else:
            fail("non-square S-box round-trip failed")
        # --- }
        subsection("Surjective: F_3^2 to F_3^1")
        # --- {
        raw_lut_psi = [[0], [1], [2], [1], [2], [0], [2], [0], [1]]
        lut_psi = [[Fp(x) for x in row] for row in raw_lut_psi]
        psi = get_sbox(lut_psi)
        if psi.get_input_size() == 2 and psi.get_output_size() == 1:
            success("psi: F_3^2 -> F_3^1 has correct dimensions (t=2, u=1)")
        else:
            fail(
                "psi dimensions wrong: t={}, u={}".format(
                    psi.get_input_size(), psi.get_output_size()
                )
            )
        # --- }
        # --- {
        if not psi.is_invertible():
            success("psi is correctly non-invertible (input_size != output_size)")
        else:
            fail("psi should not be invertible")
        # --- }
        # --- {
        if (
            fp_eq(psi[[0, 0]], [0])
            and fp_eq(psi[[1, 0]], [1])
            and fp_eq(psi[[2, 1]], [0])
        ):
            success("psi(x0,x1) = x0+x1 mod 3 evaluated correctly")
        else:
            fail("psi evaluation failed")
        # --- }
        # --- {
        psi_rt = get_sbox(psi.to_bytes())
        if isinstance(psi_rt, S_box_fp) and psi_rt == psi:
            success("non-square S-box (t=2, u=1) survives bytes round-trip")
        else:
            fail("non-square S-box (t=2, u=1) round-trip failed")
        # --- }
        subsection("Composing across different shapes")
        # --- {
        psi_phi = psi * phi
        if [list(y) for y in psi_phi] == [[0], [0], [0]]:
            success("psi * phi is the zero function (x + 2x = 0 mod 3)")
        else:
            fail(
                "psi * phi gave {}, expected all zeros".format(
                    [list(y) for y in psi_phi]
                )
            )
        # --- }
        # --- {
        u_phi = u * phi
        if u_phi.get_input_size() == 1 and u_phi.get_output_size() == 2:
            success("u * phi: F_3^1 -> F_3^2 has correct dimensions")
        else:
            fail(
                "u * phi has wrong dimensions: t={}, u={}".format(
                    u_phi.get_input_size(), u_phi.get_output_size()
                )
            )
        # --- }
        # --- {
        expected_u_phi = [list(u[list(row)]) for row in raw_lut_phi]
        if [list(y) for y in u_phi] == expected_u_phi:
            success("(u * phi)[x] == u(phi(x)) for all x in F_3")
        else:
            fail(
                "u * phi gave {}, expected {}".format(
                    [list(y) for y in u_phi], expected_u_phi
                )
            )
        # --- }
        section("Error handling")
        subsection("Incompatible sizes in arithmetic")
        # --- {
        try:
            _ = u + f
            fail("adding S-boxes with different lengths should raise")
        except Exception:
            success("u + f raises when lengths differ (9 vs 5)")
        # --- }
        # --- {
        try:
            _ = u - f
            fail("subtracting S-boxes with different lengths should raise")
        except Exception:
            success("u - f raises when lengths differ")
        # --- }
        subsection("Incompatible sizes in composition")
        # --- {
        try:
            _ = u * psi
            fail("composing u (input_size=2) with psi (output_size=1) should raise")
        except Exception:
            success("u * psi raises: output_size(psi)=1 != input_size(u)=2")
        # --- }
        # --- {
        try:
            _ = phi * u
            fail("composing phi (input_size=1) with u (output_size=2) should raise")
        except Exception:
            success("phi * u raises: output_size(u)=2 != input_size(phi)=1")
        # --- }
        subsection("Inverting a non-bijection")
        # --- {
        try:
            _ = g.inverse()
            fail("g.inverse() should raise for the non-bijective x^2 mod 5")
        except Exception:
            success("g.inverse() correctly raises for non-bijective g")
        # --- }
        # --- {
        try:
            _ = g**-1
            fail("g ** -1 should raise for non-bijective g")
        except Exception:
            success("g ** -1 correctly raises for non-bijective g")
        # --- }
        subsection("component with wrong vector length")
        # --- {
        try:
            _ = u.component([1, 0, 0])
            fail("component([1,0,0]) on output_size=2 S-box should raise")
        except Exception:
            success("component([1,0,0]) on output_size=2 S-box correctly raises")
        # --- }
        subsection("coordinate out of range")
        # --- {
        try:
            _ = u.coordinate(2)
            fail("coordinate(2) on output_size=2 S-box should raise")
        except Exception:
            success("coordinate(2) on output_size=2 S-box correctly raises")
        # --- }
        section("Mathematical properties")
        subsection("Commutativity and associativity of addition")
        # --- {
        if (u + v) == (v + u):
            success("addition is commutative: u + v == v + u")
        else:
            fail("addition is not commutative!")
        # --- }
        # --- {
        if ((u + v) + u) == (u + (v + u)):
            success("addition is associative: (u + v) + u == u + (v + u)")
        else:
            fail("addition is not associative!")
        # --- }
        # --- {
        zero_sb = u - u
        if (u + zero_sb) == u:
            success("u + zero == u: zero is the additive identity")
        else:
            fail("u + zero != u")
        # --- }
        subsection("Identity law for composition")
        # --- {
        if (u * id_sb) == u and (id_sb * u) == u:
            success("u * Id == u and Id * u == u")
        else:
            fail("identity law for composition failed")
        # --- }
        subsection("Higher power")
        # --- {
        if (u**3) == (u * u * u):
            success("u ** 3 == u * u * u")
        else:
            fail("u ** 3 != u * u * u")
        # --- }
        subsection("Inverse undoes iteration")
        # --- {
        if (u**2).inverse() == (u.inverse() ** 2):
            success("(u ** 2).inverse() == u.inverse() ** 2")
        else:
            fail("(u ** 2).inverse() != u.inverse() ** 2")
        # --- }
        subsection("Derivative in multiple directions")

        # --- {
        def ref_deriv(sbox, delta):
            Fp_s = GF(sbox.get_p())
            lut_s = [list(y) for y in sbox.get_lut()]
            sp = [list(x) for x in sbox.get_input_space()]
            powers = sbox.get_powers_in()
            in_dim = sbox.get_input_size()
            out_dim = sbox.get_output_size()
            result = []
            for i, x in enumerate(sp):
                xd = [Fp_s(x[j]) + Fp_s(delta[j]) for j in range(in_dim)]
                idx_xd = S_box_fp.vec_to_int(xd, powers)
                result.append(
                    [Fp_s(lut_s[idx_xd][j]) - Fp_s(lut_s[i][j]) for j in range(out_dim)]
                )
            return result

        # --- }
        # --- {
        delta_01 = [0, 1]
        got_01 = [list(y) for y in u.derivative(delta_01)]
        exp_01 = ref_deriv(u, delta_01)
        if got_01 == exp_01:
            success("derivative in direction (0,1) matches reference")
        else:
            fail(
                "derivative in direction (0,1): got {}, expected {}".format(
                    got_01, exp_01
                )
            )
        # --- }
        # --- {
        delta_22 = [2, 2]
        got_22 = [list(y) for y in u.derivative(delta_22)]
        exp_22 = ref_deriv(u, delta_22)
        if got_22 == exp_22:
            success("derivative in direction (2,2) matches reference")
        else:
            fail(
                "derivative in direction (2,2): got {}, expected {}".format(
                    got_22, exp_22
                )
            )
        # --- }
        subsection("Derivative of the identity gives delta")
        # --- {
        delta_12 = [1, 2]
        d_id = id_sb.derivative(delta_12)
        if [list(y) for y in d_id] == [delta_12] * 9:
            success("derivative of identity is the constant delta=(1,2) for all x")
        else:
            fail(
                "derivative of identity gave {}, expected constant {}".format(
                    [list(y) for y in d_id], [delta_12] * 9
                )
            )
        # --- }
        subsection("Linearity of the derivative in S")
        # --- {
        delta_lin = [1, 0]
        du = u.derivative(delta_lin)
        dv = v.derivative(delta_lin)
        dupv = (u + v).derivative(delta_lin)
        if dupv == (du + dv):
            success(
                "D_delta(u+v) == D_delta(u) + D_delta(v): derivative is linear in S"
            )
        else:
            fail("derivative linearity failed")
        # --- }
        subsection("Derivative of the zero function is zero")
        # --- {
        d_zero_fn = zero_sb.derivative([2, 1])
        if [list(y) for y in d_zero_fn] == [[0, 0]] * 9:
            success("derivative of the zero function is zero for any direction")
        else:
            fail(
                "derivative of zero function is not zero: {}".format(
                    [list(y) for y in d_zero_fn]
                )
            )
        # --- }
        subsection("Derivation is commutative")
        # --- {
        da = [1, 0]
        db = [0, 1]
        d_a_then_b = u.derivative(da).derivative(db)
        d_b_then_a = u.derivative(db).derivative(da)
        if d_a_then_b == d_b_then_a:
            success("second derivative is symmetric: D_a D_b(S) == D_b D_a(S)")
        else:
            fail("second derivative not symmetric!")
        # --- }
        subsection("Coordinate of derivative equals derivative of coordinate")
        # --- {
        delta_cd = [1, 0]
        deriv_coord_0 = u.derivative(delta_cd).coordinate(0)
        coord_0_deriv = u.coordinate(0).derivative(delta_cd)
        if deriv_coord_0 == coord_0_deriv:
            success(
                "D_delta(u)_0 == D_delta(u_0): derivative commutes with coordinate extraction"
            )
        else:
            fail("D_delta(u)_0 != D_delta(u_0)")
        # --- }
        subsection("Equality is reflexive, symmetric, and consistent with iteration")
        # --- {
        u3 = get_sbox(lut)
        if u == u3 and u3 == u:
            success("equality is symmetric")
        else:
            fail("equality is not symmetric")
        # --- }
        # --- {
        if u == u:
            success("equality is reflexive")
        else:
            fail("equality is not reflexive!")
        # --- }
        # --- {
        u_shifted = u + u - u  # should equal u
        if u_shifted == u:
            success("u + u - u == u: arithmetic round-trips to equality")
        else:
            fail("u + u - u != u")
        # --- }
        # --- {
        if (u + u) != u:
            success("u + u != u for non-trivial u")
        else:
            fail("u + u == u, which should not happen for non-zero u")
        # --- }
        section("DDT properties")
        # --- {
        from sboxU.statistics import (
            fp_ddt,
            fp_differential_spectrum,
            fp_differential_uniformity,
        )

        # --- }
        # --- {
        id3 = get_sbox([[Fp(0)], [Fp(1)], [Fp(2)]])
        sq3 = get_sbox([[Fp(0)], [Fp(1)], [Fp(1)]])
        # --- }
        subsection("Row sums equal input-space size")
        # --- {
        ddt_id = fp_ddt(id3)
        if all(sum(row) == 3 for row in ddt_id):
            success("every DDT row of id3 sums to 3")
        else:
            fail(
                "some DDT row of id3 has wrong sum: {}".format([sum(r) for r in ddt_id])
            )
        # --- }
        subsection("Delta=0 row")
        # --- {
        if ddt_id[0] == [3, 0, 0]:
            success("DDT row 0 of id3 is [3,0,0]")
        else:
            fail("DDT row 0 of id3 is {}, expected [3,0,0]".format(ddt_id[0]))
        # --- }
        subsection("Identity S-box rows")
        # --- {
        if ddt_id[1] == [0, 3, 0]:
            success("DDT row 1 of id3 is [0,3,0]")
        else:
            fail("DDT row 1 of id3 is {}, expected [0,3,0]".format(ddt_id[1]))
        # --- }
        # --- {
        if ddt_id[2] == [0, 0, 3]:
            success("DDT row 2 of id3 is [0,0,3]")
        else:
            fail("DDT row 2 of id3 is {}, expected [0,0,3]".format(ddt_id[2]))
        # --- }
        subsection("PN function: uniform DDT rows")
        # --- {
        ddt_sq = fp_ddt(sq3)
        if ddt_sq[1] == [1, 1, 1] and ddt_sq[2] == [1, 1, 1]:
            success("DDT of squaring map has uniform rows [1,1,1] for nonzero delta")
        else:
            fail("DDT of squaring map: row1={}, row2={}".format(ddt_sq[1], ddt_sq[2]))
        # --- }
        subsection("Differential spectrum")
        # --- {
        sp_id = fp_differential_spectrum(id3)
        from sboxU.core.spectrum.cython_functions import Spectrum

        if isinstance(sp_id, Spectrum):
            success("fp_differential_spectrum returns a Spectrum instance")
        else:
            fail(
                "fp_differential_spectrum returned {}, expected Spectrum".format(
                    type(sp_id)
                )
            )
        # --- }
        # --- {
        if sp_id[3] == 2:
            success("fp_differential_spectrum(id3)[3] == 2")
        else:
            fail("fp_differential_spectrum(id3)[3] == {}, expected 2".format(sp_id[3]))
        # --- }
        # --- {
        if sp_id[0] == 4:
            success("fp_differential_spectrum(id3)[0] == 4")
        else:
            fail("fp_differential_spectrum(id3)[0] == {}, expected 4".format(sp_id[0]))
        # --- }
        subsection("Differential uniformity")
        # --- {
        if fp_differential_uniformity(id3) == 3:
            success("differential uniformity of id3 is 3")
        else:
            fail(
                "differential uniformity of id3 is {}, expected 3".format(
                    fp_differential_uniformity(id3)
                )
            )
        # --- }
        # --- {
        if fp_differential_uniformity(sq3) == 1:
            success("differential uniformity of squaring map (PN) is 1")
        else:
            fail(
                "differential uniformity of sq3 is {}, expected 1".format(
                    fp_differential_uniformity(sq3)
                )
            )
        # --- }
    return exit_code()


if __name__ == "__main__":
    sys.exit(main_test())
