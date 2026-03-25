from sage.all import *
from sboxU import *
from sboxU.display import *
from sage.crypto.sboxes import sboxes


# n = 8
# s = random_permutation_S_box(8)
# s = F2_trans(s[0], bit_length=8) * s
# print(differential_spectrum(s))
# print(thickness_spectrum(s))
# le_repr = le_class_representative(s)
# print(le_repr)
# print(differential_spectrum(le_repr))
# print(linear_equivalence_permutations(s, le_repr))
# print(linear_equivalence_permutations(s, random_permutation_S_box(8)))


# !SECTION! Specific Tests


def test_ccz_exploration():
    with Experiment("testing CCZ-exploration"):
        n = 6
        cube = monomial(3, GF(2**n))
        print(thickness_spectrum(cube), cube.get_input_length(), cube.get_output_length())
        ws = get_WalshZeroesSpaces_quadratic_apn(cube)
        print(ws.thickness_spectrum())
        for L in ws.get_mappings():
            g = ccz_equivalent_function(cube, L)
            ws_prime = ws.image_by(L.transpose().inverse())
            print("")
            pprint(thickness_spectrum(g))
            pprint(ws_prime.thickness_spectrum())

    
def test_affine_equivalence():
    with Experiment("Testing the affine equivalence of Midori_Sb0"):
        
        Sb0=get_sbox(sboxes["Midori_Sb0"])
        for c in ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']:
            S_test = get_sbox(sboxes["Optimal_S"+c])
            res = affine_equivalence(Sb0, S_test)
            if len(res)>0:
                print("Sb0 is affine equivalent to Optimal_S"+c,'\n')
                subsection("Trying to reconstruct Sb0 from the affine_equivalence result")
                T_out = F2_trans(res[3],bit_length=4)
                L_out = get_sbox(res[2])
                T_in = F2_trans(res[1],bit_length=4)
                L_in = get_sbox(res[0])
                if T_out*L_out*S_test*L_in*T_in == Sb0:
                    print("Success")
                else :
                    print("Failure")
            else:
                print("not equivalent")
        

def test_EA_mappings(verbose=False):
    with Experiment("Testing EA-mappings"):
        n = 6
        field = GF(2**n)
        alpha = field.gen()
        X = PolynomialRing(field, "X").gen()
        A = Blm(alpha * X)
        B = Blm(X**2)
        C = Blm(X**4 + alpha**3*X)
        
        cube = X**3
        
        if verbose:
            print(A)
            print("")
            print(B)
            print("")
            print(C)
            print("")
        L = EA_mapping(A, B, C)
        s0 = ccz_equivalent_function(cube, L)
        s1 = get_sbox(B)*cube*get_sbox(A) + get_sbox(C)
        if verbose:
            print(L)
            print(s0)
            print(s1)
        if s0 == s1:
            print("Success")
        else:
            print("Failure")


def test_ccz_equivalence_to_permutation(verbose=False):
    with Experiment("Testing EA-mappings"):
        n = 6
        field = GF(2**n)
        alpha = field.gen()
        X = PolynomialRing(field, "X").gen()
        kim_mapping = get_sbox(X**3 + alpha*X**24 + X**10)
        # ws = get_WalshZeroesSpaces_quadratic_apn(kim_mapping)
        for p in enumerate_permutations_in_ccz_class(kim_mapping):
            print(is_permutation(p))
            pprint(thickness_spectrum(p), degree_spectrum(p))
        

# !SECTION! Main function 

if __name__ == "__main__":
    test_ccz_exploration()
    test_affine_equivalence()
    test_EA_mappings()
    test_ccz_equivalence_to_permutation()
