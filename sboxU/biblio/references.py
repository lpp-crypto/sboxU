import re
from pybtex.database import parse_file
from os.path import dirname

"""The path to the bibtex file containing all entries. """
bib_file_path = dirname(__file__) + "/biblio.bib"


SboxU_BIBLIOGRAPHY = None

def open_bibliography():
    global SboxU_BIBLIOGRAPHY
    if SboxU_BIBLIOGRAPHY == None:
        SboxU_BIBLIOGRAPHY = parse_file(bib_file_path)
    else:
        pass
    

def cite(key):
    """Parses bib file shipped with sboxU and returns the bibtex entry
    corresponding to the given key.

    Throws an error if there is no entry with the given key.

    """
    open_bibliography()
    target = key
    if key[0] == "[":
        target = target[1:]
    if key[-1] == "]":
        target = target[:-1]
    if target in SboxU_BIBLIOGRAPHY.entries.keys():
        return SboxU_BIBLIOGRAPHY.entries[target]
    else:
        return "{} not found in bibliography".format(target)


        
def who_to_cite(sboxU_tool):
    """Returns a list of strings, each being a bibtex entry citing a paper on which the object `sboxU_tool` is based.

    It works by finding the line "refs: [key0], [key1]..." in the docstring, and parsing it to extract all the keys, i.e. all the substrings of the form "[conf-abbrev:namesYY]", where "conf-abbrev" is a string containing only letters, "names" is also a string containing only letters, and YY is a sequence of two digits.

    Args:
        sboxU_tool: a function or a module of sboxU. Anything with a docstring in fact.

    !WARNING! We need to do something to be able to handle papers which should have identical keys
    """
    return re.findall(r"\[.+:.+[0-9]\]", sboxU_tool.__doc__ )


def cite_as(key, output_format):
    entry = cite(key)
    if isinstance(entry, (str)):
        return '**WARNING** bibtex entry "{}" was not found'.format(key)
    else:
        if output_format in ["text", "yaml"]:
            return entry.to_string(output_format)
        elif output_format == "rst":
            return format_to_rst(entry, key)
            
        

def format_to_rst(entry, key):
    content = entry.fields
    content["key"] = key
    for person_type in [ "author", "editor" ]:
        if person_type in entry.persons.keys():
            content[person_type] = ""
            for p in entry.persons[person_type]:
                for x in p.first_names:
                    content[person_type] += x + " "
                for x in p.last_names:
                    content[person_type] += x + " "
                content[person_type] = content[person_type][:-1] + ", "
            content[person_type] = content[person_type][:-2]
    return "\\[{key}\\]\n    {author} ({year}). *{title}*. `\\[link\\] <{url}>`_".format(**entry.fields)


def gen_full_biblio_rst(path):
    open_bibliography()
    with open(path, "w") as f:
        f.write(".. This file is auto-generated, modifications are pointless!\n")
        f.write("Bibliography\n=======\n\n")
        for key in sorted(SboxU_BIBLIOGRAPHY.entries.keys()):
            f.write("\n")
            f.write(cite_as(key, "rst"))
        f.write("\n")
