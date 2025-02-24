#!/usr/bin/env sage
# Time-stamp: <2025-02-24 14:27:39>


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

# simply download the `py` at
#
# https://github.com/lpp-crypto/levain/tree/master
#
# and rename it `levain`
from levain import *

from sage.crypto.sboxes import sboxes
from sboxU.known_functions import *


HASH_ALGORITHM = hashlib.sha256
ID_BITLENGTH   = 128


# !SECTION! S-box encoding as bytes
# =================================


def get_block_lengths(s):
    """Return the number of bits `n` in the input and `m` in the
    output of `s`

    """
    # finding n
    n = 1
    while (1 << n) < len(s):
        n += 1
    if 2**n != len(s):
        raise Exception("wrong S-box length")
    else:
        # finding m
        mask = 1
        m = 1
        for x in s:
            while (x != (x & mask)):
                mask = (mask << 1) | 1
                m += 1
        return n, m


def pack_to_bytes(s, m):
    """Packs the content of s into a sequence of 8-bit blocks, i.e. a
    `bytearray`. Assumes that all elements of `s` are strictly smaller
    than 2**m.

    """
    if m <= 4:
        result = [0]*floor(len(s) / 2)
        for i in range(0, len(s), 2):
            result[i >> 1] = (s[i] << 4) | s[i+1]
        if (len(s) % 2) == 1: # handling the case of an odd length
            result.append(s[-1])
        return bytearray(result)
    elif m <= 8:
        return bytearray(s)
    else:
        byte_length = ceil(m / 8)
        result = [0] * len(s) * byte_length
        for i in range(0, len(s)):
            x = s[i]
            for j in range(0, byte_length):
                result[i * byte_length + j] = x & 0xFF
                x = x >> 8
        return bytearray(result)

    
def encode_lut(s, m):
    """Returns an array of bytes encoding the full lookup table `s`.

    Since tinySQL doesn't support arrays, we instead store them as
    BLOBs.

    """
    return  bytes([m]) + pack_to_bytes(s, m)
    

def decode_lut(l):
    """Returns a list of integers corresponding to the lut encoded by
    the bytearray l.

    """
    m = l[0]
    b = [int(x) for x in l[1:]]
    if m <= 4:
        result = [0] * len(b)
        for x in range(0, len(result), 2):
            y = b[(x >> 1)]
            result[x]   = y >> 4
            result[x+1] = y & 0xF
    else:
        block = ceil(m / 8)
        result = [0] * int(len(b) / block)
        for x in range(0, len(result)):
            result[x] = sum(Integer(b[x*block + j]) << (8*j) for j in range(0, block))
    return result
    


# !SECTION! S-box identifiers
# ===========================


def hash_as_integer(b):
    digest = HASH_ALGORITHM(b).digest()
    result = 0
    for x in digest[0:ceil(ID_BITLENGTH/8)]:
        result = 256*result + x
    return result


def apn_identifier(lut):
    try:
        o = ortho_derivative(lut)
        spectras = "{} || {}".format(
            pretty_spectrum(differential_spectrum(o)),
            pretty_spectrum(walsh_spectrum(o), absolute=True)
        )
        return hash_as_integer(spectras.encode("UTF-8"))
    except:
        n, m = get_block_lengths(lut)
        return hash_as_integer(pack_to_bytes(lut, m))





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
    


# !SECTION! Wrappers
# ==================
                 

            
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
                    if content > 0:  # case of positive query on integer
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
        inserted_list = [entry[k] for k in sorted(self.row_structure.keys())]
        try:
            self.cursor.execute(self.function_insertion_query, tuple(inserted_list))
        except:
            raise Exception("Insertion failed for \n {}\n".format(entry))


    # handling the `with` syntax
        
    def __enter__(self):
        self.connection = sqlite3.connect(self.db_file)
        self.cursor = self.connection.cursor()
        return self

    def __exit__(self, *args):
        self.connection.commit()
        self.connection.close()


        
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

class WalshZeroesSpaces:
    def __init__(self, blob=None, lut=None):
        if blob == None and lut == None:
            raise Exception("need at least a blob representation or the lut of the original function!")
        elif lut == None:
            self.init_from_blob(blob)
        else:
            self.init_from_lut(lut)
        self.thickness_spectrum = thickness_spectrum([], spaces=self.spaces, N=self.n)
        self.Ls = [get_generating_matrix(V, 2*self.n).transpose()
                   for V in self.spaces]
    

    def init_from_lut(self, lut):
        self.n, self.m = get_block_lengths(lut)
        self.spaces = get_lat_zeroes_spaces(lut)

        
    def init_from_blob(self, b):
        self.n = int(b[0])
        self.m = int(b[1])
        b = b[2:]
        block = ceil((self.m+self.n) / 8)
        bases = [0] * int(len(b) / block)
        for x in range(0, len(bases)):
            bases[x] = sum(Integer(b[x*block + j]) << (8*j) for j in range(0, block))
        self.spaces = [bases[i:i+self.n] for i in range(0, len(bases), self.n)]


    def to_blob(self):
        all_spaces = []
        for v in self.spaces:
            all_spaces += v
        return bytes([self.n, self.m]) + pack_to_bytes(all_spaces, self.n + self.m)


    def reduce(self, automorphisms):
        relevant = [True] * len(self.Ls)
        for a in ELEMENTS_OF(range(0, len(self.Ls)), "admissible mappings"):
            La_inv = self.Ls[a].inverse()
            for b in range(a+1, len(self.Ls)):
                if relevant[b]:
                    for u in automorphisms:
                        if is_EA(self.Ls[b] * u * La_inv):
                            relevant[b] = False
                            break
        self.spaces = [self.spaces[a]
                       for a in range(0, len(self.spaces)) if relevant[a]]
        self.Ls     = [self.Ls[a]
                       for a in range(0, len(self.Ls)) if relevant[a]]


    


    # !TODO! an __rmult__ function to apply a FastLinearMapping to it
    
    def __str__(self):
        return "Walsh spaces with thicknesses {}".format(pretty_spectrum(self.thickness_spectrum))
    
        
    def __len__(self):
        return len(self.spaces)


    def __iter__(self):
        for v in self.spaces:
            yield v
    
    

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
                "walsh_spaces": "BLOB"
            }
        )

        
    def insert_function_from_lut(self, lut, bibliography):
        n, m = get_block_lengths(lut)
        encoded = encode_lut(lut, n)
        # linear
        walsh_spec = walsh_spectrum(lut)
        lin = 0
        for k in walsh_spec.keys():
            lin = max(lin, abs(k))
        # degree
        deg = max([a.degree() for a in algebraic_normal_form(lut)])
        # thickness and spaces
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
            "walsh_spaces" : spaces.to_blob()
        }
        self.insert_function(to_insert)

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
            # post-processing
        entry["lut"] = decode_lut(entry["lut"])
        entry["walsh_spaces"] = WalshZeroesSpaces(blob=entry["walsh_spaces"])
        return entry

    
        
def test_APNFunctions():
    with LogBook("Testing APNFunctions"):

        
        SECTION("filling the database")
        with APNFunctions() as db:
            db.create()
            from sboxU.known_functions import sixBitAPN as reservoir_6
            from sboxU.known_functions import sevenBitAPN as reservoir_7
            for s in ELEMENTS_OF(reservoir_6.all(), "6 bits"):
                db.insert_function_from_lut(s, "Banff")
            for s in ELEMENTS_OF(reservoir_7.all(), "7 bits"):
                db.insert_function_from_lut(s, "Dunno")

    
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



    
# !SECTION! Testing
# =================

def is_EA(M):
    half_lines = int(M.nrows() / 2)
    half_cols  = int(M.ncols() / 2)
    for k, l in itertools.product(range(0, half_lines),
                                  range(half_cols, M.ncols())):
        if M[k][l] != 0:
            return False
    return True
                

def is_identity(M):
    for i in range(0, M.nrows()):
        for j in range(0, M.ncols()):
            if (i == j) and (M[i][j] != 1):
                return False
            elif (i != j) and (M[i][j] != 0):
                return False
    return True



if __name__ == "__main__":
    # test_LiteratureSBoxes()
    
    test_APNFunctions()
