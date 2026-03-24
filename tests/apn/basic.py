from sage.all import *
from sboxU import *



with Experiment("Testing APN functions-related functions"):

    section("Initialization")
    n = 7
    g = GF(2**n)
    cube = monomial(3, g)
    print("n=", n)
    print(cube)

    section("Testing mugshot for quadratic APN functions")
    pprint(differential_spectrum(cube))
    pprint(ortho_derivative(cube))
    pprint(sigma_multiplicities(cube, 4))
    pprint(thickness_spectrum(cube))
    
    print(apn_ea_mugshot(cube))

    section("Testing CCZ-class reduction")
    
    reprs = enumerate_ea_classes_apn_quadratic(cube)
    print("EA classes in CCZ: ", len(reprs))
    for s in reprs:
        pprint(degree_spectrum(s))
        print(s)
        pprint("checking: ", degree_spectrum(ccz_equivalent_quadratic_function(s)))

    section("Testing DB")
    
    with APNFunctions("./apnDB.db") as db:
        s = db[0]["sbox"]
        pprint(s)
        pprint(differential_spectrum(s))
        pprint(thickness_spectrum(s))
        pprint(degree_spectrum(s))
        pprint(absolute_walsh_spectrum(s))

