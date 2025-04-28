#!/usr/bin/env sage
# Time-stamp: <2025-04-28 17:25:53>


# /!\ You are not really supposed to look at this file: highly
# experimental stuff is happening! S-boxes database *are* coming to
# sboxU, but you shouldn't build upon or rely in the content of this
# file in any way: it is not stable at all, and my undergo
# API-breaking changes without any notice.

from sage.all import *
from collections import defaultdict

from sboxU import *
from switching import *

import base64
import hashlib
import sqlite3

# simply download the `py` folder at
#
# https://github.com/lpp-crypto/levain/tree/master
#
# and rename it `levain`
from levain import *

from sage.crypto.sboxes import sboxes
from sboxU.known_functions import *



# !SECTION! S-box identifiers
# ===========================



def mugshot(s,
            degree_spec=None,
            walsh_spec=None,
            differential_spec=None,
            thickness_spec=None):
    """Returns a string such that, if this string is different for two
    functions, then they cannot be EA-equivalent.

    It can be seen as a "sketch" of the function: if two functions
    have the same mugshot, then they might be EA-equivalent. If not,
    they definitely aren't. There are mutliple ways to achieve this
    depending on the specifics of the function:
    
    - for a quadratic APN function, we use the differential and Walsh
      spectra of their ortho-derivative

    - otherwise, we use the concatenation of the differential,
      absolute Walsh, thickness, and degree spectra.

    """
    if degree_spec == None:
        degree_spec = degree_spectrum(s)
    if differential_spec == None:
        differential_spec = differential_spectrum(s)
    if max(degree_spec.keys()) == 2 and max(differential_spec.keys()) == 2:
        # case of a quadratic APN
        o = ortho_derivative(s)
        return "Quad APN; w{} d{}".format(
            pretty_spectrum(walsh_spectrum(o)), # note that we don't take absolute values
            pretty_spectrum(differential_spectrum(o))
        )
    else:
        if walsh_spec == None:
            walsh_spec = walsh_spectrum(s)
        if thickness_spec == None:
            thickness_spec = thickness_spectrum(s)
        result = "w{} d{} deg{} thk{}".format(
            pretty_spectrum(walsh_spec, absolute=True),
            pretty_spectrum(differential_spec),
            pretty_spectrum(degree_spec),
            pretty_spectrum(thickness_spec),
        )
        if max(differential_spec.keys()) == 2:
            result += "sig{}".format(pretty_spectrum(sigma_multiplicities(s, 4)))
        return result



# !SECTION! Data base generation
# ==============================




def db_sboxes_setup():
    """Creates the "literature_sboxes.db" file and stores the AES,
    Stribog, and PRESENT S-boxes in it.

    """
    print("* Writing")
    print("** initializing table")
    connection = sqlite3.connect("literature_sboxes.db")
    cursor = connection.cursor()
    creation_query = "CREATE TABLE IF NOT EXISTS functions ("
    for column in sorted(ROW_STRUCTURE_SBOX.keys()):
        creation_query += "{} {},".format(column, ROW_STRUCTURE_SBOX[column])
    creation_query = creation_query[:-1] + ")"
    cursor.execute(creation_query)
        
    print("** inserting S-boxes")
    insertion_query = "INSERT INTO functions VALUES ("
    insertion_query += "?, " * len(ROW_STRUCTURE_SBOX)
    insertion_query = insertion_query[:-2] + ")"

    aes_sbox = {
        "lut" : encode_lut(list(sboxes["AES"]), 8),
        "n" : 8,
        "m" : 8,
        "bibliography" : "DaeRij99",
        "differential_uniformity" : 4,
        "linearity" : 32,
        "degree" : 7
    }

    stribog_sbox = {
        "lut" : encode_lut(list(sboxes["Stribog"]), 8),
        "n" : 8,
        "m" : 8,
        "bibliography" : "Russia",
        "differential_uniformity" : 8,
        "linearity" : 56,
        "degree" : 7
    }
    

    present_sbox = {
        "lut" : encode_lut(list(sboxes["PRESENT"]), 4),
        "n" : 4,
        "m" : 4,
        "bibliography" : "presentPaper",
        "differential_uniformity" : 4,
        "linearity" : 8,
        "degree" : 3
    }

    twine_sbox = {
        "lut" : encode_lut(list(sboxes["TWINE"]), 4),
        "n" : 4,
        "m" : 4,
        "bibliography" : "twinePaper",
        "differential_uniformity" : 4,
        "linearity" : 8,
        "degree" : 3
    }
    
    for entry in [aes_sbox, stribog_sbox, twine_sbox, present_sbox]:
        inserted_list = [entry[k] for k in sorted(ROW_STRUCTURE_SBOX.keys())]
        cursor.execute(insertion_query, tuple(inserted_list))
    connection.commit()
    connection.close()
    


# !SECTION! Main Data base
# ========================
                 

class FunctionDB:
    """This idea of this class is to factor away the interaction with
    any database of S-boxes (i.e. with both the database of S-boxes
    from the literature and with the APN function database).

    In particular, it handles the generation of the SELECT queries,
    and the (admitedly small) boilerplate needed by the `with` syntax.

    """

    # !TODO! handle the bibliography
    # !TODO! how to store spectras in a way that allows meaningful queries? 
    
    def __init__(self, db_file, row_structure):
        self.db_file = db_file
        self.row_structure = row_structure
        self.row_structure["id"] = "INTEGER"
        self.functions_table = "functions"
        self.bibliography_table = "bibliography"
        # preparing queries
        self.function_insertion_query = "INSERT INTO {} VALUES ({} ?)".format(
            self.functions_table,
            "?, " * (len(row_structure) - 1) # the last question mark
                                             # is already in the
                                             # string above
        )

        
    def create(self):
        creation_query = "CREATE TABLE IF NOT EXISTS {} (".format(self.functions_table)
        for column in sorted(self.row_structure.keys()):
            creation_query += "{} {},".format(column, self.row_structure[column])
        creation_query = creation_query[:-1] + ")"
        self.cursor.execute(creation_query)
        self.number_of_functions = 0
        
        
    # handling searches

    
    def parse_function_from_row(self, row):
        raise Exception("virtual method that shouldn't be called")

    
    def query_functions(self, query_description, max_entries=-1):
        where_clause = ""
        for constraint in query_description.keys():
            if constraint not in self.row_structure.keys():
                raise Exception("unrecognized parameter in query : {} (={})".format(
                    constraint,
                    query_description[constraint]
                    ))
            else:
                content = query_description[constraint]
                # case of an integer equality or difference
                if isinstance(content, (int, Integer)):
                    if content >= 0:  # case of positive query on integer
                        where_clause += " ({}=={:d}) AND".format(constraint, content)
                    else:
                        where_clause += " ({}!={:d}) AND".format(constraint, - content)
                # case of a range
                elif isinstance(content, (type(range(0,1)))):
                    where_clause += " ({}>={:d} AND {}<{:d}) AND".format(
                        constraint,
                        content.start,
                        constraint,
                        content.stop
                    )
                # case of a string
                elif isinstance(content, (str)):
                    if "%" in content:
                        # case of search
                        where_clause += " ({} like '{}') AND".format(
                            constraint,
                            content
                        )
                    else:
                        # case of an exact match
                        where_clause += " ({}=='{}') AND".format(
                            constraint,
                            content
                        )
        where_clause = where_clause[:-4] # removing the final "AND"
        self.cursor.execute("SELECT * FROM {} WHERE {}".format(
            self.functions_table,
            where_clause
        ))
        result = []
        for row in self.cursor.fetchall():
            result.append(self.parse_function_from_row(row))
            if max_entries > 0 and len(result) == max_entries:
                return result
        return result
        
    
    # handling insertions

    def insert_function(self, entry):
        entry["id"] = self.number_of_functions
        inserted_list = [entry[k] for k in sorted(self.row_structure.keys())]
        try:
            self.cursor.execute(self.function_insertion_query, tuple(inserted_list))
            self.number_of_functions += 1
            return entry["id"]
        except:
            raise Exception("Insertion failed for \n {}\n".format(entry))


    def __len__(self):
        return self.number_of_functions
    

    # handling the `with` syntax
        
    def __enter__(self):
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
        # if the file exists, we initialize the length; otherwise we don't
        try:
            self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
            self.number_of_functions = self.cursor.fetchall()[0][0]
        except:
            self.number_of_functions = 0
        return self

    def __exit__(self, *args):
        self.connection.commit()
        self.connection.close()


# !SECTION! Wrappers
# ==================
        
# !SUBSECTION! S-boxes from the cryptography literature

class LiteratureSBoxes(FunctionDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "literature_sboxes.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains just a bit of
    logic on top to handle the computation of all the S-box properties
    we are interested in.

    # !TODO! add representative of linear equivalence classes; and the
    # !logic to use it

    """
    
    def __init__(self):
        super().__init__(
            "literature_sboxes.db",
            {
                "lut" : "BLOB",
                "n" : "INTEGER",
                "m" : "INTEGER",
                "differential_uniformity" : "INTEGER",
                "linearity" : "INTEGER",
                "degree" : "INTEGER",
                "bibliography" : "TEXT",
                "cipher": "TEXT"
            }
        )

        
    def insert_function_from_lut(self, lut, cipher, bibliography):
        n, m = get_block_lengths(lut)
        encoded = encode_lut(lut, n)
        # differential
        diff_spec = differential_spectrum(lut)
        diff_unif = max(diff_spec.keys())
        # linear
        walsh_spec = walsh_spectrum(lut)
        lin = 0
        for k in walsh_spec.keys():
            lin = max(lin, abs(k))
        # degree
        deg = max([a.degree() for a in algebraic_normal_form(lut)])
        # building the query
        to_insert = {
            "lut" : encoded,
            "n" : n,
            "m" : m,
            "bibliography" : bibliography,
            "differential_uniformity" : diff_unif,
            "linearity" : lin,
            "degree" : deg,
            "cipher" : cipher
        }
        self.insert_function(to_insert)

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
            # post-processing
        entry["lut"] = decode_lut(entry["lut"])
        return entry

        
        
def test_LiteratureSBoxes():
    with LogBook("Testing LiteratureSboxes"):

        SECTION("filling the database")
        with LiteratureSBoxes() as db:
            db.create()
            db.insert_function_from_lut(list(sboxes["AES"]),    "AES",    "DaeRij99")
            db.insert_function_from_lut(list(sboxes["Stribog"]),"Stribog","Russia")
            db.insert_function_from_lut(list(sboxes["PRESENT"]),"PRESENT","presentPaper")
            db.insert_function_from_lut(list(sboxes["TWINE"]),  "TWINE",  "twinePaper")
    
        SECTION("retrieving content")
        with LiteratureSBoxes() as db:
            SUBSECTION("8-bit S-boxes")
            for entry in db.query_functions({"n" : 8}):
                print(entry)
            SUBSECTION("Differentially-4 S-boxes")
            for entry in db.query_functions({"differential_uniformity" : 4}):
                print(entry)
        
    


# !SUBSECTION! APN vectorial Boolean functions

    

class APNFunctions(FunctionDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "apn_functions.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains additional
    logic to handle the specifics of APN functions, and in particular
    of their CCZ-equivalence class. Here, "functions" should be
    thought of much more as "extended affine equivalence class
    representative" rather than function.
    
    This class provides another table containing the bases of all the
    spaces of dimension n contained with the Walsh zeroes of a
    function. In order to both save space and store the structure of a
    CCZ-equivalence, each APN function is stored along with the
    identifier of the Walsh spaces of its CCZ-equivalence class, and
    the FastLinearMapping that must be applied to it to obtain its own
    Walsh Zeroes.

    """
    
    def __init__(self):
        # !IDEA! have a max_degree and a min_degree?
        super().__init__(
            "apn_functions.db",
            {
                "lut" : "BLOB", # the lookup table of the representative
                "n" : "INTEGER",
                "m" : "INTEGER",
                "linearity" : "INTEGER",
                "thickness" : "INTEGER",
                "degree" : "INTEGER",
                "bibliography" : "TEXT",
                "ccz_id": "INTEGER",
                "walsh_L": "BLOB",
                "mugshot" : "TEXT"
            }
        )
        self.spaces_table = "spaces"
        self.spaces_insertion_query = "INSERT INTO {} VALUES (?, ?)".format(
            self.spaces_table,
        )

    def __enter__(self):
        super().__enter__()
        try:
            self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.spaces_table))
            self.number_of_ccz_classes = self.cursor.fetchall()[0][0]
        except:
            self.number_of_ccz_classes = 0
        return self

    def __str__(self):
        return "APN function DB containing {} EA-classes from {} CCZ-classes".format(
            self.number_of_functions,
            self.number_of_ccz_classes
        )

            
    def create(self):
        # creating the functions table
        super().create()
        # spaces table name
        spaces_creation_query  = "CREATE TABLE IF NOT EXISTS {}".format(self.spaces_table)
        # spaces table format
        spaces_creation_query +=  "(id INTEGER, bases BLOB)"
        # creating the table
        self.cursor.execute(spaces_creation_query)
        self.number_of_ccz_classes = 0


    def insert_function_from_lut(self, lut, bibliography):
        differential_spec = differential_spectrum(lut)
        if max(differential_spec.keys()) != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(lut))
        n, m = get_block_lengths(lut)
        encoded = encode_lut(lut, n)
        # linear
        walsh_spec = walsh_spectrum(lut)
        lin = 0
        for k in walsh_spec.keys():
            lin = max(lin, abs(k))
        # degree
        degree_spec = degree_spectrum(lut)
        deg = max(degree_spec.keys())
        # thickness and spaces
        spaces = WalshZeroesSpaces(lut=lut)
        thk_spec = spaces.thickness_spectrum()
        thk = max(thk_spec.keys())
        # inserting the spaces
        self.cursor.execute(self.spaces_insertion_query, (
            self.number_of_ccz_classes,
            spaces.to_blob()
        ))
        self.number_of_ccz_classes += 1
        # inserting the function 
        to_insert = {
            "lut" : encoded,
            "n" : n,
            "m" : m,
            "bibliography" : bibliography,
            "linearity" : lin,
            "degree" : deg,
            "thickness" : thk,
            "ccz_id" : self.number_of_ccz_classes,
            "walsh_L" : encode_lut(list(range(0, n+m)), n+m), # the identity
            "mugshot" : mugshot(lut,
                                walsh_spec=walsh_spec,
                                degree_spec=degree_spec,
                                differential_spec=differential_spec,
                                thickness_spec=thk_spec)
        }
        return self.insert_function(to_insert)

    
    def fetch_WalshZeroesSpaces(self, index):
        self.cursor.execute("SELECT * FROM {} WHERE id={:d}".format(
            self.spaces_table,
            index
        ))
        for row in self.cursor.fetchall():
            W = WalshZeroesSpaces(blob=row[1])
            return W
        raise Exception("no WalshZeroesSpaces with index={}".format(index))


    def insert_full_ccz_equivalence_class(self, lut, bibliography):
        differential_spec = differential_spectrum(lut)
        if max(differential_spec.keys()) != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(lut))
        n, m = get_block_lengths(lut)
        encoded = encode_lut(lut, n)
        # linear
        walsh_spec = walsh_spectrum(lut)
        lin = 0
        for k in walsh_spec.keys():
            lin = max(lin, abs(k))
        # finding WalshZeroesSpaces
        spaces = WalshZeroesSpaces(lut=lut)
        if algebraic_degree(lut) == 2: # if the function is quadratic,
                                       # we compute automorphisms
                                       # inserting the spaces
            quadratic_case = True
            autom = [FastLinearMapping(L)
                     for L in graph_automorphisms_of_apn_quadratic(lut)]
            Ls = spaces.reduced_linear_mappings(autom)
        else:
            quadratic_case = False
            Ls = spaces.Ls
        self.cursor.execute(self.spaces_insertion_query, (
            self.number_of_ccz_classes,
            spaces.to_blob()
        ))
        # inserting all the functions
        for L in Ls:
            new_lut = apply_mapping_to_graph(lut, L)
            if quadratic:
                worth_adding = True
            else:
                worth_adding = (self.is_present(new_lut, test_ccz=True)[0] == "absent")
            if worth_adding:
                L_inv_T = L.transpose().inverse()
                new_spaces = L_inv_T * spaces
                new_thk_spec = new_spaces.thickness_spectrum()
                # t1 = pretty_spectrum(thickness_spectrum(new_lut))
                # t2 = pretty_spectrum(new_thk_spec)
                # if t1 != t2:
                #     FAIL("well there's your problem: {} != {}".format(t1, t2))
                new_degree_spec =degree_spectrum(new_lut)
                to_insert = {
                    "lut" : encode_lut(new_lut, n),
                    "n" : n,
                    "m" : m,
                    "bibliography" : bibliography,
                    "linearity" : lin,
                    "degree" : max(new_degree_spec.keys()),
                    "thickness" : max(new_thk_spec.keys()),
                    "ccz_id" : self.number_of_ccz_classes,
                    "walsh_L" : encode_lut(L_inv_T.masks, n+m),
                    "mugshot" : mugshot(lut,
                                        walsh_spec=walsh_spec,
                                        degree_spec=new_degree_spec,
                                        differential_spec=differential_spec,
                                        thickness_spec=new_thk_spec)
                }
                self.insert_function(to_insert)
        self.number_of_ccz_classes += 1
        return self.number_of_ccz_classes
    

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        entry["lut"] = decode_lut(entry["lut"])
        W = self.fetch_WalshZeroesSpaces(entry["ccz_id"])
        L = FastLinearMapping(decode_lut(entry["walsh_L"]))
        entry["walsh_spaces"] = L * W
        return entry

    
    def is_present(self,
                   s,
                   degree_spec=None,
                   walsh_spec=None,
                   differential_spec=None,
                   thickness_spec=None,
                   test_ccz=True):
        """Returns one of three things:
        
        - `["present", index]` if there is already a function
          extended-affine equivalent to `s` in the data-base where
          `index` is its `id`,
        
        - `["absent"]` if there isn't, and
        
        - [`"maybe"`, index]` if at least one function with a similar
          mugshot is present, that function having an `id` equal to
          `index`.

        If `test_ccz` is set to `True` (the default case), then a
        CCZ-equivalence test is performed using SAGE's implementation
        of code equivalence to get a definitive answer. `"maybe"` can
        only be returned if `test_ccz` is set to `False`.

        """
        n, m = get_block_lengths(s)
        encoded = encode_lut(s, n)
        mug = mugshot(s,
                      degree_spec=degree_spec,
                      walsh_spec=walsh_spec,
                      differential_spec=differential_spec,
                      thickness_spec=thickness_spec)
        candidates = self.query_functions({"mugshot" : mug})
        if len(candidates) == 0:
            return ["absent"]
        elif test_ccz:
            # checking all the functions with a similar mugshot
            for entry in candidates:
                if entry["lut"] == encoded:
                    return ["present", entry["id"]]
                else:
                    #print("manually testing CCZ equivalence")
                    if are_ccz_equivalent(entry["lut"], s):
                        return ["present", entry["id"]]
            # if we checked for CCZ-equivalence and didn't find any
            # CCZ-equivalent function, then `s` is new...
            return ["absent"]
        else:
            # ... but if we didn't check, then we can't be sure.
            return ["maybe", candidates[0]["id"]]

        
    
# !SECTION! Testing
# =================
    

# !SUBSECTION! Testing APN function database


def test_APNFunctions():
    with LogBook("Testing APNFunctions"):

        SECTION("filling the database")
        with APNFunctions() as db:
            db.create()
            from sboxU.known_functions import sixBitAPN as reservoir_6
            from sboxU.known_functions import sevenBitAPN as reservoir_7
            for s in ELEMENTS_OF(reservoir_6.all_quadratics()[4:], "6 bits"):
                db.insert_function_from_lut(s, "Banff")
    
        SECTION("retrieving content")
        with APNFunctions() as db:
            SUBSECTION("high linearity")
            for entry in db.query_functions({"linearity": range(20, 34)}):
                print(entry, desc="*l")
                v = WalshZeroesSpaces(lut=entry["lut"])
            SUBSECTION("high thickness")
            for entry in db.query_functions({"thickness": range(0, 2)}):
                print(entry["walsh_spaces"], desc="l*")
                w = WalshZeroesSpaces(lut=entry["lut"])
                print("Checking if Walsh zeroes storing works ({})".format(w))
                b = w.to_blob()
                w_prime = WalshZeroesSpaces(blob=b)
                if (w.spaces == w_prime.spaces):
                    SUCCESS("encoding works")
                else:
                    FAIL("mismatch between direct computation and decoding")


def test_APNFunctions_CCZ_insert():
    with LogBook("Testing APNFunctions CCZ insertion"):

        SECTION("filling the database")
        with APNFunctions() as db:
            db.create()
            print(db)
            from sboxU.known_functions import sixBitAPN as reservoir_6
            for s in ELEMENTS_OF(reservoir_6.all_quadratics(), "6 bits"):
                db.insert_full_ccz_equivalence_class(s, "CCZ")
    
        SECTION("retrieving content")
        with APNFunctions() as db:
            print(db)
            SUBSECTION("high linearity")
            for entry in db.query_functions({"linearity": range(16, 34)}):
                print(entry["id"], desc="*l")
            SUBSECTION("high thickness")
            for entry in db.query_functions({"thickness": range(4, 6)}):
                print(entry["ccz_id"], desc="l*")
                u = entry["walsh_spaces"]
                w = WalshZeroesSpaces(lut=entry["lut"])
                print("Checking if Walsh zeroes storing works ({})".format(u))
                b = w.to_blob()
                w_prime = WalshZeroesSpaces(blob=b)
                if (w.spaces == w_prime.spaces):
                    SUCCESS("encoding works")
                else:
                    FAIL("mismatch between direct computation and decoding")
                if (w.spaces == u.spaces):
                    SUCCESS("storing works")
                else:
                    FAIL("mismatch between direct computation and recovery from DB")
                    print("{}  {}".format(pretty_spectrum(w.thickness_spectrum()),
                                          pretty_spectrum(u.thickness_spectrum())))



# !SECTION! Initializing the data-bases
# =====================================


def fill_6bit_APNFunctions():
    with LogBook("Filling the APN database"):
        with APNFunctions() as db:
            SECTION("Initialization")
            db.create()
            from sboxU.known_functions import sixBitAPN as reservoir
            SECTION("Case of quadratic functions")
            for index, s in enumerate(reservoir.all_quadratics()):
                SUBSECTION("function n°{}".format(index))
                db.insert_full_ccz_equivalence_class(s, "Banff")
            SECTION("Case of the non-CCZ quadratic function")
            unique_non_quadratic = [0, 0, 0, 8, 0, 26, 40, 58, 0, 33, 10, 35, 12, 55, 46, 29, 0, 11, 12, 15, 4, 21, 32, 57, 20, 62, 18, 48, 28, 44, 50, 10, 0, 6, 18, 28, 10, 22, 48, 36, 8, 47, 16, 63, 14, 51, 62, 11, 5, 24, 27, 14, 11, 12, 61, 50, 25, 37, 13, 57, 27, 61, 39, 9]
            db.insert_full_ccz_equivalence_class(unique_non_quadratic, "non-quadratic")
                    

# !SUBSECTION! Main function


def print_func(s):
    print("lut={}\ndeg={}\nthk={}\nsig={}".format(
        s,
        pretty_spectrum(degree_spectrum(s)),
        pretty_spectrum(thickness_spectrum(s)),
        pretty_spectrum(sigma_multiplicities(s, 4))
    ))
    

if __name__ == "__main__":
    # test_LiteratureSBoxes()
    
    # test_APNFunctions()
    
    # test_APNFunctions_CCZ_insert()

    # fill_6bit_APNFunctions()

    # with LogBook("Inspecting new CCZ-class"):
    #     with APNFunctions() as db:
    #         SECTION("Grabbing new functions")
    #         funcs = db.query_functions({"id" : range(780, 900)})
    #         for entry in funcs:
    #             s = entry["lut"]
    #             SECTION("function with id={}".format(entry["id"]))
    #             print("lut={}\ndeg={}\nthk={}".format(
    #                 s,
    #                 pretty_spectrum(degree_spectrum(s)),
    #                 pretty_spectrum(thickness_spectrum(s))
    #             ))
    
    with LogBook("Exploration 6-bit"):
        j = None
        with APNFunctions() as db:
            SECTION("Grabbing original APN functions")
            print(db)
            funcs = db.query_functions({"n": 6, "degree": range(3, 8)})
            shuffle(funcs)
            for entry in funcs:
                SECTION("APN function with id={}".format(entry["id"]), timed=True)
                s = entry["lut"]
                print_func(s)
                SUBSECTION("computing switching neighbours", timed=True)
                candidates = apn_non_affine_switching_class(s)
                print("# candidates found: {}".format(len(candidates)))
                index = 0
                shuffle(candidates)
                seen = defaultdict(int)
                for swi in ELEMENTS_OF(candidates[:50], "switching neighbours"):
                    is_known = db.is_present(swi, test_ccz=True)
                    if is_known[0] == "absent":
                        SUBSECTION("neighbour n°{}".format(index))
                        swi_id = db.insert_function_from_lut(swi, "switching")
                        SUCCESS("found!")
                        print_func(swi)
                        swi_entry = db.query_functions({"id": swi_id})[0]
                        w = swi_entry["walsh_spaces"]
                        for L in w.Ls:
                            f = apply_mapping_to_graph(swi, L)
                            if db.is_present(f, test_ccz=True)[0] == "absent":
                                j = db.insert_function_from_lut(f)
                                print("added new function with id={}".format(j))
                                print_func(f)
                            else:
                                FAIL("somehow we got a known CCZ equivalence class from an unknown function")
                    else:
                        print("known ({})".format(is_known[1]), desc="n")
                        seen[is_known[1]] += 1
                        if seen[is_known[1]] > 20:
                            break
                    index += 1


    # with APNFunctions() as db:
    #     entries = db.query_functions({"degree": 4})
    #     for entry in entries:
    #         print("deg={} ; thk={}".format(
    #             pretty_spectrum(degree_spectrum(entry["lut"])),
    #             pretty_spectrum(thickness_spectrum(entry["lut"]))
    #         ))
