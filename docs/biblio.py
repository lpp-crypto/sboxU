#!/usr/bin/env sage

from sboxUv2 import *

bibliography_path = "./source/bibliography.rst"
print("writing to " + bibliography_path + "...")
gen_full_biblio_rst(bibliography_path)
print("[DONE]")
