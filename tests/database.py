from sage.all import *
from sboxUv2 import *

from sage.crypto.sboxes import sboxes

def test_LiteratureSBoxes():
    with Experiment("Testing LiteratureSboxes"):

        section("filling the database")
        
        with LiteratureSBoxes() as db:
            db.create()
            print("DB created")
            
            db.insert_function_from_lut(
                sboxes["AES"],
                "AES",
                "DaeRij99",
                usage="AES,SNOW"
            )
            db.insert_function_from_lut(
                sboxes["Stribog"],
                "Stribog",
                "GOST",
                usage="Stribog,Kuznyechik,Kuznechik"
            )
            db.insert_function_from_lut(
                sboxes["PRESENT"],
                "PRESENT",
                "presentPaper"
            )
            db.insert_function_from_lut(
                sboxes["TWINE"],
                "TWINE",
                "twinePaper"
            )
    
        section("retrieving content")
        
        with LiteratureSBoxes() as db:
            subsection("8-bit S-boxes")
            for entry in db.query_functions({"n" : 8}):
                print(entry)
            subsection("Differentially-4 S-boxes")
            for entry in db.query_functions({"differential_uniformity" : 4}):
                print(entry)
            subsection("Grabbing the PRESENT S-box")
            for entry in db.query_functions({"lut" : sboxes["PRESENT"]}):
                print(entry)
            subsection("Grabbing the Kuznyechik S-box")
            for entry in db.query_functions({"usage" : "%Kuznechik%"}):
                print(entry)
        

if __name__ == "__main__":
    test_LiteratureSBoxes()
