from sboxUv2 import *
from sys import argv

DEFAULT_APN_DB_NAME = "apnDB_{}.db"

def generate_apn_ea_classes_database(
        n,
        db_path=None
):
    """!TODO! write docstring
    """
    with Experiment("Generating database for n={}".format(n)):
        section("Grabbing class representatives")

        subsection("reading")
        
        # handling inputs
        if n == 6:
            from reprs6 import ccz_class_representatives
        elif n == 7:
            from reprs7 import ccz_class_representatives
        if db_path == None:
            db_path = DEFAULT_APN_DB_NAME.format(n)
        # generating the DB
        print("path = ", db_path)
        print("read {} functions".format(len(ccz_class_representatives)))
        subsection("splitting into quadratic and non-quadratic")
        quads = []
        non_quads = []
        for s in ccz_class_representatives:
            if algebraic_degree(s) == 2:
                quads.append(s)
            else:
                non_quads.append(s)
    

        with APNFunctions(db_path) as db:
    
            section("filling the databse with CCZ-quadratic EA classes")
            
            for index, s in enumerate(quads):
                inserted = db.insert_full_ccz_equivalence_class(
                    s,
                    "6-bit"
                )
                print("{:5d}) CCZ class has {} functions".format(index, len(inserted)))
    
            section("filling the database with non-CCZ-quadratic EA classes")

            for s in non_quads:
                pprint(thickness_spectrum(s))
                pprint(degree_spectrum(s))
                inserted = db.insert_full_ccz_equivalence_class(
                    s,
                    "{}-bit".format(n)
                )
                print("CCZ class with {} functions added".format(len(inserted)))
                
            print("generation finished")
            print("{} EA classes generated".format(len(db)))
            print("{} CCZ classes found".format(db.number_of_ccz_classes))


if __name__ == "__main__":
    n = int(argv[1])
    if n not in [6,7]:
        print("n must be 6 or 7")
    else:
        generate_apn_ea_classes_database(n, "apn-{}.db".format(n))
