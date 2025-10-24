from sboxUv2 import *
from pathlib import Path
import argparse


DEFAULT_APN_DB_NAME = "apnDB.db"


def generate_apn_ea_classes_database(
        ccz_class_representatives,
        db_path
):
    """!TODO! write docstring
    """
    with Experiment("Generating database"):
        
        section("Grabbing class representatives")

        subsection("reading")
        
        # generating the DB
        print("path = ", db_path)
        if Path(db_path).is_file():
            print("[WARNING] file already exists, we are deleting it")
            Path(db_path).unlink()
        else:
            print("DB will be stored in {}".format(db_path))
        
        print("using {} class representatives".format(len(ccz_class_representatives)))
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
                    "--"
                )
                print("CCZ class with {} functions added".format(len(inserted)))
                
            print("generation finished")
            print("{} EA classes generated".format(len(db)))
            print("{} CCZ classes found".format(db.number_of_ccz_classes))

            section("How to use it?")

            print("\n1. Put the file {} in the folder with your script".format(db_path))
            print("\n2. make sure you import sboxUv2 in your sage script")
            print("\n3. use the following syntax to check if F is new:")
            print("with APNFunctions('{}') as db:".format(db_path))
            print("    print('Yay!' if db.is_new(F) else ':((')")
            print("\n4. use the following to print the Walsh spectrum of all degree 3 functions:")
            print("with APNFunctions('{}') as db:".format(db_path))
            print("    for entry in db.query_functions({'degree' : 3}):")
            print("        s = entry['sbox']")
            print("        pprint(walsh_spectrum(s))")
            print("\n\n")

def process_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path",
                        nargs="?",
                        type=str,
                        default=DEFAULT_APN_DB_NAME,
                        help="The path to the file in which to write the database")
    parser.add_argument("-n", 
                        type=int,
                        help="The value of the bit length of the functions to process")
    return parser.parse_args()

    
def main_cli():
    args = process_arguments()
    n = args.n
    path = args.path
    print(n, path)
    if n == 6:
        from .reprs6 import ccz_class_representatives
    elif n == 7:
        from .reprs7 import ccz_class_representatives
    if n not in [6,7]:
        print("n must be 6 or 7")
    else:
        generate_apn_ea_classes_database(
            ccz_class_representatives,
            path
        )


if __name__ == "__main__":
    main_cli()
