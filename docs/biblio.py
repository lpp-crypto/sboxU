#!/usr/bin/env sage

from sboxUv2 import *
from sys import argv

bibliography_rst  = "./source/bibliography.rst"
bibliography_html = "./build/html/bibliography.html"

def gen_bibliography_dot_rst():
    print("writing to " + bibliography_rst + "...")
    gen_full_biblio_rst(bibliography_rst)
    print("[DONE]")


def add_ids():
    new_content = []
    with open(bibliography_html, "r") as f:
        for line in f.readlines():
            result = re.findall(r">\[.+:.+[0-9]\]", line)
            # !CONTINUE! replace the string by one that includes an "id"


# !TODO! write a function that parses sboxUv2.[a-z]*.html files and replaces bib keys with links 

if __name__ == "__main__":
    if argv[1] == "gen":
        gen_bibliography_dot_rst()
