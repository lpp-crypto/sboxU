from sage.all import *
from sboxUv2 import *
from sage.crypto.sboxes import sboxes


n = 8
s = random_permutation_S_box(8)
s = F2_trans(s[0], bit_length=8) * s
print(differential_spectrum(s))
print(thickness_spectrum(s))
le_repr = le_class_representative(s)
print(le_repr)
print(differential_spectrum(le_repr))
print(linear_equivalence_permutations(s, le_repr))
print(linear_equivalence_permutations(s, random_permutation_S_box(8)))


# !SECTION! Testing CCZ exploration

n = 7
cube = monomial(3, GF(2**n))
print(thickness_spectrum(cube), cube.get_input_length(), cube.get_output_length())
s_s = enumerate_ea_classes(cube)
print("tot: ", len(s_s))
for s in s_s:
    print(s)#, thickness_spectrum(s))

# !SECTION! Testing affine equivalence on 4 bits

print('\n',"----Testing the affine equivalence of Midori_Sb0----",'\n')
Sb0=Sb(sboxes["Midori_Sb0"])
for c in ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']:
    S_test=Sb(sboxes["Optimal_S"+c])
    res=affine_equivalence(Sb0,S_test)
    if len(res)>0:
        print("Sb0 is affine equivalent to Optimal_S"+c,'\n')
        print("----Trying to reconstruct Sb0 from the affine_equivalence result----")
        if F2_trans(res[3],bit_length=4)*res[2].get_S_box()*S_test*res[0].get_S_box()*F2_trans(res[1],bit_length=4)==Sb0:
            print("Success")
        else :
            print("Failure")

