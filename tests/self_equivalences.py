from sboxU import *
from sage.crypto.sboxes import sboxes

if __name__ == "__main__":
    with Experiment("Testing self_equivalences functions"):
        section("Testing if the mappings are correct")
        subsection("For Midori_Sb0's Sbox (4-bit)")

        S=get_sbox(sboxes["Midori_Sb0"])

        for elt in self_affine_equivalent_mappings(S):
            print(elt[0].get_S_box()*S*elt[1].get_S_box()==S)

        # subsection("For Scream's Sbox (8-bit)")
        # print("The computation might take some time ...")
        # S=get_sbox(sboxes["Scream"])

        # for elt in self_affine_equivalent_mappings(S):
        #     print(elt[0].get_S_box()*S*elt[1].get_S_box()==S)






