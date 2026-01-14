from sboxUv2 import *
from pathlib import Path
import argparse
import hashlib


########################################################
#### !SECTION! Code For Generic EA Classes Database ####
########################################################

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
    
            section("filling the database with CCZ-quadratic EA classes")
            
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


################################################################
#### !SECTION! Code For CCZ only quadratic compact database ####
################################################################


def are_ea_equivalent_from_vq(f,g):

    # !! TODO !!
    # Add base cases

    ws_f = get_WalshZeroesSpaces(f)
    mappings_f = list(ws_f.get_mappings())


    ws_g = get_WalshZeroesSpaces(g)
    mappings_g = list(ws_g.get_mappings())

    

    quad_index_f = 0
    for i in range(len(mappings_f)):
        f_eq = ccz_equivalent_function(f, mappings_f[i])
        if algebraic_degree(f_eq) == 2:
             quad_index_f = i
             break

    quad_index_g = 0
    for i in range(len(mappings_g)):
        g_eq = ccz_equivalent_function(g, mappings_g[i])
        if algebraic_degree(g_eq) == 2:
             quad_index_g = i
             break
    
    
    # Compute of the ccz-quadratic aut(q)
    q = ccz_equivalent_quadratic_function(f)
    aut_q = automorphisms_from_ortho_derivative(q)
    
    # Compute the EA mapping
    q_prime = ccz_equivalent_quadratic_function(g)
    ea_mappings = ea_mappings_from_ortho_derivative(q_prime,q)
    
    if len(ea_mappings) != 0 : 
        ea_q_q_prime = ea_mappings[0]
    else:
         return(False)

    # To check
    #print(q_prime == ccz_equivalent_function(q,ea_q_q_prime))

    ws_1 =  ws_f.image_by(mappings_f[quad_index_f].inverse().transpose())
    ws_1 =  ws_1.image_by(mappings_f[quad_index_f].inverse().transpose())

    ws_2 =  ws_g.image_by(mappings_g[quad_index_g].inverse().transpose())
    ws_2 =  ws_2.image_by(mappings_g[quad_index_g].inverse().transpose())
    ws_2 =  ws_2.image_by(ea_q_q_prime.transpose())

    #print(ws_1.get_bases()[quad_index_f])
    #print(ws_2.get_bases()[quad_index_g])


    for L in aut_q:
        ws_temp = ws_1.image_by(L.transpose().inverse())
        if ws_2.get_bases()[quad_index_g]  == ws_temp.get_bases()[quad_index_f]:
            return(True)

    return(False)

DEFAULT_APN_DB_NAME = "apnDB.db"

class APNQuadraticFunctions_ccz_only(FunctionsDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "apn_functions.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains additional
    logic to handle the specifics of APN functions, and in particular
    of their CCZ-equivalence class. Here, "functions" should be
    thought of much more as "extended affine equivalence class
    representative" rather than function.
    
    This class is a compact version of APNFunctions to use when there are
    too many functions when we include EA-classes. 

    """
    
    def __init__(self, db_file):
        # !IDEA! have a max_degree and a min_degree?
        super().__init__(
            db_file,
            {
                "qcr" : "BLOB",
                "mugshot" : "BLOB"
            }
        )
        if self.new_db:
            self.number_of_ccz_classes = 0
        else:
            try:
                self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
            except:
                self.number_of_ccz_classes = 0
                
    
    def __str__(self):
        return "APN function DB containing {} CCZ-classes".format(
            self.number_of_functions

        )


    def insert_quadratic_ccz_representative(self, s):
        """
        In this version, we only insert a ccz representative
        """
        
        sb = Sb(s)
        differential_spec = differential_spectrum(sb)
        if differential_spec.maximum() != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(sb))

        if algebraic_degree(sb) == 2: 
            mug = apn_ea_mugshot(sb)
            # We hash the mugshot for memory concerns
            h = hashlib.sha256()
            h.update(mug)
            mug = h.digest()

            worth_adding = self.is_new(sb,mug)
        else:
            s_quad = ccz_equivalent_quadratic_function(sb.lut())
            if s_quad != []:
                sb = s_quad
                mug = apn_ea_mugshot(sb)
                # We hash the mugshot for memory concerns
                h = hashlib.sha256()
                h.update(mug)
                mug = h.digest()
                worth_adding = self.is_new(sb,mug)
            else:
                raise Exception("Trying to add a non ccz_quadratic APN function to the database: {}".format(sb))



        if worth_adding:
            to_insert = {
                "qcr" : Sb(quadratic_compact_representation(sb.lut())).to_bytes(),
                "mugshot" : mug
            }
            inserted_id = self.insert_function(to_insert)
            self.number_of_ccz_classes += 1
            return inserted_id
        return None
    

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        entry["sbox"] = Sb(quadratic_sbox_from_compact_representation(entry["qcr"],8,8))
        return entry

    
    def is_new(self,
                   s,
                   mug=None
        ):
        """
        # !TODO! docstring 
        """
        sb = Sb(s)
        
        if mug == None:
            # Computing the mugshot of s depending on its degree
            if degree_spec == None:
                degree_spec = degree_spectrum(sb)
            if degree_spec.maximum() == 2:
                mug = apn_ea_mugshot(sb)
            else:
                # If s is ccz quadratic, we work with a quadratic representative
                s_quad = ccz_equivalent_quadratic_function(sb)
                if s_quad == []:
                    # !! TODO !!
                    # Decide if we put more
                    mug = absolute_walsh_spectrum(sb)
                else:
                    mug = apn_ea_mugshot(Sb(s_quad))


            # We hash the mugshot for memory concerns
            h = hashlib.sha256()
            h.update(mug)
            mug = h.digest()


        query = {"mugshot" : mug}
        candidates = self.query_functions(query)
        if candidates == []:
            #print("- New Mugshot")
            return True
        else:
            for entry in candidates:
                print("Same Mugshot as Function {}".format(entry["id"]))
                print(entry["mugshot"])
                
                # Trivial Case
                if quadratic_compact_representation(sb.lut()) == entry["qcr"]:
                    return False
                
                # Testing EA-equivalence
                qcr = entry["qcr"]
                lut = quadratic_sbox_from_compact_representation(Sb(qcr).lut(),sb.get_input_length(),sb.get_output_length())
                b = are_ea_equivalent_from_vq(sb.lut(), lut)
                if b:
                    print("Function EA equivalent to  id = {} in the database".format(entry["id"]))
                    return False

            return True


def generate_apn_ccz_classes_database(
        ccz_class_representatives,
        db_path
):
    """!TODO! write docstring
    """
    with Experiment("Generating database"):
        
        section("Grabbing class representatives")

        subsection("reading")
        
        # Generating the DB
        print("path = ", db_path)
        if Path(db_path).is_file():
            print("[WARNING] file already exists, we are deleting it")
            Path(db_path).unlink()
        else:
            print("DB will be stored in {}".format(db_path))
        
        print("using {} class representatives".format(len(ccz_class_representatives)))
        with APNQuadraticFunctions_ccz_only(db_path) as db:
    
            section("filling the database with CCZ-quadratic classes")
            
            for index, s in enumerate(ccz_class_representatives):
                
                #print("Trying {}th Function in Representative List".format(index))
                
                if algebraic_degree(s) == 2:
                    inserted = db.insert_quadratic_ccz_representative(s)
                    if inserted == None :
                        print("- {}th Function is not New".format(index))
                        print()
                    #else:
                    #    print("- {}th Function is New".format(index))
                    #    print()
                else:
                    print("Non Quadratic Function, discarded")
                    print()

            print("generation finished")
            print("{} CCZ classes found".format(db.number_of_ccz_classes))

 
def main_cli():
    args = process_arguments()
    n = args.n
    path = args.path
    mode = 'ea'
    #mode  = 'ccz_compact'
    print(n, path,mode)
    if n == 6:
        from .reprs6 import ccz_class_representatives
    elif n == 7:
        from .reprs7 import ccz_class_representatives
    elif n == 8:
        from .reprs8 import ccz_class_representatives
    if n not in [6,7,8]:
        print("n must be 6,7 or 8")
    else:

        if mode == 'ccz_compact':

            generate_apn_ccz_classes_database(
                ccz_class_representatives,
                path
                    )
        else:

            generate_apn_ea_classes_database(
                ccz_class_representatives,
                path
            )


if __name__ == "__main__":
    main_cli()
