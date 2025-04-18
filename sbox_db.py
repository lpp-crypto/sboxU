#!/usr/bin/env sage
# Time-stamp: <2025-04-18 17:34:13>


# /!\ You are not really supposed to look at this file: highly
# experimental stuff is happening! S-boxes database *are* coming to
# sboxU, but you shouldn't build upon or rely in the content of this
# file in any way: it is not stable at all, and my undergo
# API-breaking changes without any notice.

from sage.all import *
from collections import defaultdict

from sboxU import *

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
        return "w{} d{} deg{} thk{}".format(
            pretty_spectrum(walsh_spec, absolute=True),
            pretty_spectrum(differential_spec),
            pretty_spectrum(degree_spec),
            pretty_spectrum(thickness_spec),
        )
    



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
        
        
    # handling searches

    
    def parse_function_from_row(self, row):
        raise Exception("virtual method that shouldn't be called")

    
    def query_functions(self, query_description):
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
        return [self.parse_function_from_row(row) 
                for row in self.cursor.fetchall()]

    
    # handling insertions

    def insert_function(self, entry):
        entry["id"] = len(self)
        inserted_list = [entry[k] for k in sorted(self.row_structure.keys())]
        try:
            self.cursor.execute(self.function_insertion_query, tuple(inserted_list))
            return entry["id"]
        except:
            raise Exception("Insertion failed for \n {}\n".format(entry))


    def __len__(self):
        self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
        return self.cursor.fetchall()[0][0]
    

    # handling the `with` syntax
        
    def __enter__(self):
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
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

    It builds upon the `FunctionDB` class, and contains just a bit of
    logic on top at this stage.

    """
    
    def __init__(self):
        # !IDEA! have a max_degree and a min_degree?
        super().__init__(
            "apn_functions.db",
            {
                "lut" : "BLOB",
                "n" : "INTEGER",
                "m" : "INTEGER",
                "linearity" : "INTEGER",
                "thickness" : "INTEGER",
                "degree" : "INTEGER",
                "bibliography" : "TEXT",
                "walsh_spaces": "BLOB",
                "mugshot" : "TEXT"
            }
        )

    def insert_function_from_lut(self, lut, bibliography, spaces=None):
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
        if spaces == None:
            spaces = WalshZeroesSpaces(lut=lut)
        thk_spec = spaces.thickness_spectrum
        thk = max(thk_spec.keys())
        # building the query
        to_insert = {
            "lut" : encoded,
            "n" : n,
            "m" : m,
            "bibliography" : bibliography,
            "linearity" : lin,
            "degree" : deg,
            "thickness" : thk,
            "walsh_spaces" : spaces.to_blob(),
            "mugshot" : mugshot(lut,
                                walsh_spec=walsh_spec,
                                degree_spec=degree_spec,
                                differential_spec=differential_spec,
                                thickness_spec=thk_spec)
        }
        return self.insert_function(to_insert)
    

    def insert_ccz_equivalent_function(self, index, L):
        entry = self.query_functions({"id": index})[0]
        s = entry["lut"]
        if not isinstance(L, FastLinearMapping):
            L = FastLinearMapping(L)
        lut = apply_mapping_to_graph(s, L)
        w = entry["walsh_spaces"]
        L_T = L.transpose()
        L_inv = L.inverse()
        for i in range(0, len(w)):
            w.Ls[i] = L * w.Ls[i]
            img = [L_T(x) for x in linear_span(w.spaces[i])]
            img.sort()
            w.spaces[i] = extract_basis(img, entry["n"]+entry["m"])
            
        return self.insert_function_from_lut(lut, entry["bibliography"], spaces=w)

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
            # post-processing
        entry["lut"] = decode_lut(entry["lut"])
        entry["walsh_spaces"] = WalshZeroesSpaces(blob=entry["walsh_spaces"])
        return entry

    
    def is_present(self,
                   s,
                   degree_spec=None,
                   walsh_spec=None,
                   differential_spec=None,
                   thickness_spec=None):
        """Returns one of three things:
        
        - `["present", index]` if there is already a function
          extended-affine equivalent to `s` in the data-base where
          `index` is its `id`,
        
        - `["absent"]` if there isn't, and
        
        - [`"maybe"`, index]` if at least one function with a similar
          mugshot is present, that function having an `id` equal to
          `index`.
        
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
        else:
            for entry in candidates:
                if entry["lut"] == encoded:
                    return ["present", entry["id"]]
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
            for entry in db.query_functions({"thickness": range(4, 20)}):
                print(entry["walsh_spaces"], desc="l*")
                print("Checking if blob encoding of Walsh zeroes works")
                w = WalshZeroesSpaces(lut=entry["lut"])
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
            from sboxU.known_functions import sixBitAPN as reservoir_6
            for index, s in enumerate(reservoir_6.all_quadratics()):
                SUBSECTION("function n°{}".format(index))
                print("first {}".format(s))
                first_id = db.insert_function_from_lut(s, "Banff")
                w = WalshZeroesSpaces(lut=s)
                db.insert_ccz_equivalent_function(first_id, w.Ls[3])
                f = apply_mapping_to_graph(s, w.Ls[3])
                db.insert_function_from_lut(f, "placeholder")
                print("number of functions: {}".format(len(db)))
        SECTION("retrieving content")
        with APNFunctions() as db:
            for entry in db.query_functions({"id": range(0, 100)}):
                print(entry)
                w = WalshZeroesSpaces(lut=entry["lut"])
                print(w)



# !SECTION! Initializing the data-bases 


def fill_6bit_APNFunctions():
    with LogBook("Filling the APN database"):
        with APNFunctions() as db:
            SECTION("Initialization")
            db.create()
            from sboxU.known_functions import sixBitAPN as reservoir
            SECTION("Case of quadratic functions")
            for index, s in enumerate(reservoir.all_quadratics()):
                SUBSECTION("function n°{}".format(index))
                w = WalshZeroesSpaces(lut=s)
                autom = [FastLinearMapping(L)
                         for L in graph_automorphisms_of_apn_quadratic(s)]
                first_id = db.insert_function_from_lut(s, "Banff", spaces=w)
                w.reduce(autom)
                for k, L in enumerate(w.Ls):
                    f = apply_mapping_to_graph(s, L)
                    if max(degree_spectrum(f).keys()) == 2:
                        print("ignoring: quadratic", desc="n")
                    else:
                        j = db.insert_ccz_equivalent_function(first_id, L)
                        print("inserted at 'id'={}".format(j), desc="n")
                            
                


# !SUBSECTION! Main function

if __name__ == "__main__":
    # test_LiteratureSBoxes()
    
    # test_APNFunctions_CCZ_insert()

    # fill_6bit_APNFunctions()

    unique_non_quadratic = [0, 0, 0, 8, 0, 26, 40, 58, 0, 33, 10, 35, 12, 55, 46, 29, 0, 11, 12, 15, 4, 21, 32, 57, 20, 62, 18, 48, 28, 44, 50, 10, 0, 6, 18, 28, 10, 22, 48, 36, 8, 47, 16, 63, 14, 51, 62, 11, 5, 24, 27, 14, 11, 12, 61, 50, 25, 37, 13, 57, 27, 61, 39, 9]

    w = WalshZeroesSpaces(lut=unique_non_quadratic)
    by_mug = defaultdict(list)
    dif = differential_spectrum(unique_non_quadratic)
    wal = walsh_spectrum(unique_non_quadratic)
    for index, L in enumerate(w.Ls):
        s = apply_mapping_to_graph(unique_non_quadratic, L)
        L_T = L.transpose()
        w_s = []
        for i in range(0, len(w)):
            img = [L_T(x) for x in linear_span(w.spaces[i])]
            img.sort()
            w_s.append(extract_basis(img, 12))
        print("\n", index)
        print(s)
        deg = degree_spectrum(s)
        thk = thickness_spectrum(s, spaces=w_s)
        print("lut: {}\ndeg: {}\nthk: {}".format(
            s,
            pretty_spectrum(deg),
            pretty_spectrum(thk),
        ))
        mug = mugshot(s,
                      degree_spec      = deg,
                      walsh_spec       = wal,
                      differential_spec= dif,
                      thickness_spec   = thk)
        by_mug[mug].append(s)
    for c in sorted(by_mug):
        print("\nsize = {}\nmug  = {}".format(
            len(by_mug[c]),
            c
        ))
        if len(by_mug[c]) > 1:
            fs = by_mug[c]
            kept = [True for i in range(0, len(fs))]
            for i in range(0, len(fs)-1):
                if kept[i]:
                    for j in range(i+1, len(fs)):
                        if kept[j]:
                            if are_ccz_equivalent(fs[i], fs[j]):
                                kept[j] = False
            reduced_fs = [fs[i] for i in range(0, len(fs)) if kept[i]]
            print("started from {}, filtered into {}".format(len(fs), len(reduced_fs)))
            for f in reduced_fs:
                print(f)

    # !CONTINUE! clean up the handling of the non-quadratic one and move it to the fill_6bit function 
    
    # with LogBook("CCZ-class exploration 6-bit non-quadratic"):
    #     SECTION("generating targets")
    #     funcs = known_functions.sixBitAPN.all_non_quadratics()[:]
    #     SECTION("Looping through {} functions".format(len(funcs)))
    #     for index, s in enumerate(funcs):
    #         SUBSECTION("Function n° {}".format(index))
    #         print("lut: {}".format(s))
    #         # o = ortho_derivative(s)
    #         # print("o-dif: {}\no-lin: {}".format(
    #         #     pretty_spectrum(differential_spectrum(o)),
    #         #     pretty_spectrum(walsh_spectrum(o), absolute=True)
    #         # ))
    #         if algebraic_degree(s) == 2:
    #             SUBSUBSECTION("generating automorphisms", timed=True)
    #             autom = [FastLinearMapping(L)
    #                      for L in graph_automorphisms_of_apn_quadratic(s)]
    #             print("{} automorphisms found".format(len(autom)))
    #         else:
    #             autom = None
    #         SUBSUBSECTION("Generating EA representatives", timed=True)
            
    #         reprs = enumerate_ea_classes(s, automorphisms=autom)
    #         print("number of EA-classes: {}".format(len(reprs)))
    #         SUBSUBSECTION("Listing representatives ({} found)".format(len(reprs)))
    #         for index, f in enumerate(reprs):
    #             # PARAGRAPH("EA-class {}".format(index))
    #             d = pretty_spectrum(degree_spectrum(f))
    #             # print("lut: {}\ndeg: {}\nthk: {}".format(
    #             #     f,
    #             #     d,
    #             #     pretty_spectrum(thickness_spectrum(f))
    #             # ))
    #             if d == "{2: 63}":
    #                 FAIL("CCZ-quadratic")
    #                 break
    #             elif index == len(reprs) - 1:
    #                 SUCCESS("NOT CCZ-quadratic")
