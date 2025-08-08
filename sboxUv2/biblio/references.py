import re
from pybtex.database import parse_file
from os.path import dirname

"""The path to the bibtex file containing all entries. """
bib_file_path = dirname(__file__) + "/biblio.bib"


def cite(key):
    """Parses bib file shipped with sboxU and returns the bibtex entry
    corresponding to the given key.

    Throws an error if there is no entry with the given key.

    """
    print(bib_file_path)
    library = parse_file(bib_file_path)
    return library.entries[key].to_string("bibtex")

        
def who_to_cite(sboxU_tool):
    """Returns a list of strings, each being a bibtex entry citing a paper on which the object `sboxU_tool` is based.

    It works by finding the line "refs: [key0], [key1]..." in the docstring, parsing it to extract the key, and then calling the function `cite` on each.

    Args:
        - sboxU_tool: a function or a module of sboxU. Anything with a docstring in fact.
    """
    # !TODO! write the `who_to_cite` function 
    result = []
    for line in sboxU_tool.__doc__:
        citation_line = re.search(r".*Refs:(.*)$", str)
        print(citation_line)
