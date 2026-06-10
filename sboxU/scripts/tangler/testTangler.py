from sys import argv
import re
import os

from sboxU.biblio import *

SHEBANG = "#!/usr/bin/env python"
REFERENCE_TITLE = "References"
PREAMBLE_TITLE  = "Preamble"
BLOCK_DELIMITER = "```"
NOT_A_HEADER = -1
DEFAULT_PREAMBLE = ["import sys",
                    "from sage.all import *",
                    "from sboxU import *"]


# !SECTION! The TestTangler class

def header_structure(line : str) -> int:
    """Computes the depth of the header corresponding to the given line, and returns its depth and the header itself.

    For example, returns (2, "some heading") on the input "## some heading"

    Returns If the line is a header, then (depth, "header"), otherwise (NOT_A_HEADER, "") if the line is not a header
    
    """
    if not line.startswith("#"):
        return (NOT_A_HEADER(), "")
    else:
        for i in range(0, len(line)):
            if line[i] != "#":
                return (i, line[i:].strip())
        raise Exception("something went wrong: end of line somehow reached in `header_depth`")
            

class TestTangler:
    """This class is used to parse a mardown file containing SAGE code with some minimal structure requirements, and to turn it into an sboxU-style `Experiment`.
    
    """

    PREAMBLE  = "# Tangled from {}\n\n{}\n\nwith Experiment('{}'):\n"
    INDENT    = "    "
    BLOCK_BEGIN = "# --- { \n"
    BLOCK_END   = "# --- } \n"
    FOOTER= "\nif __name__ == '__main__':    sys.exit(main_test())\n"

    
    def __init__(self, md_file_path : str, verbose : bool=False):
        self.original_file = md_file_path
        self.tangled_file  = md_file_path.replace(".md", ".py")
        self.verbose = verbose
        self.experiment_has_started = False
        self.experiment_title = "Experiment"
        self.current_block = None
        self.original = None
        self.tangled  = None
        self.references = []
        self.preamble = DEFAULT_PREAMBLE
        
        
    def finalize_current_block(self) -> None:
        if self.experiment_has_started:
            local_indent = 2*self.INDENT
            self.tangled.write(local_indent + self.BLOCK_BEGIN)
            for line in self.current_block:
                self.tangled.write(local_indent + line)
            self.tangled.write(local_indent + self.BLOCK_END)
            if self.verbose:
                print("tangled {} lines".format(len(self.current_block)))
        else:
            self.preamble += self.current_block
        self.current_block = None    

            
    def start_experiment(self) -> None:
        print(SHEBANG)
        for line in self.preamble:
            self.tangled.write(line + "\n")
        self.tangled.write("\n\ndef main_test():\n")
        self.tangled.write(self.INDENT + "with Experiment('{}'):\n".format(
            self.experiment_title
        ))
        self.experiment_has_started = True
        if self.verbose:
            print("Experiment: ", title)
        
            
    def process_section(self, line : str) -> None:
        depth, title = header_structure(line)
        if depth == NOT_A_HEADER:
            return
        elif depth == 1:
            self.experiment_title = title
            if self.verbose:
                print("experiment title found: ", self.experiment_title)
        elif depth == 2:
            if PREAMBLE_TITLE not in title:
                if self.experiment_has_started == False:
                    self.start_experiment()
                self.tangled.write(2*self.INDENT + "section('{}')\n".format(title))
        elif depth == 3:
            if PREAMBLE_TITLE not in title:
                self.tangled.write(2*self.INDENT + "subsection('{}')\n".format(title))
        else: # all other sections
            if PREAMBLE_TITLE not in title:
                self.tangled.write(2*self.INDENT + "print('{}')\n".format(title))
        
        
    def absorb_line(self, line : str) -> None:
        if line.startswith(BLOCK_DELIMITER):
            if self.current_block != None:
                self.finalize_current_block()
                if self.verbose:
                    print("--- }")
            else:
                self.current_block = []
                if self.verbose:
                    print("--- {")
        elif isinstance(self.current_block, (list)):
            self.current_block.append(line)
            if self.verbose:
                print(line, len(self.current_block))
        elif line.startswith("#"):
            self.process_section(line)
        else:
            self.maybe_add_reference(line)

            
    def maybe_add_reference(self, line : str) -> None:
        for hit in re.findall(r"\[\^[A-Za-z]+-[A-Za-z\+]+[0-9]*\]", line):
            ref = hit.replace("^", "")
            if ref not in self.references:
                self.references.append(ref)
            print(self.references[-1])
        

    def process(self):
        if self.verbose:
            print("Tangling from {}".format(self.original_file))
        with open(self.original_file, "r") as original:
            self.original = original
            with open(self.tangled_file, "w") as tangled:
                self.tangled = tangled
                for line in original.readlines():
                    if len(line) > 1:
                        self.absorb_line(line)
                        if line.rstrip().endswith(REFERENCE_TITLE):
                            break
                self.tangled.write(self.INDENT + "return exit_code()\n\n")
                self.tangled.write(self.FOOTER)
        if len(self.references) > 0:
            if self.verbose:
                print("handling refs")
            self.add_references()

            
    def add_references(self) -> None:
        tmp_file = self.original_file + ".tmp"
        with open(self.original_file, "r") as original:
            with open(tmp_file, "w") as updated:
                for line in original.readlines():
                    if line.rstrip().endswith(REFERENCE_TITLE):
                        break
                    else:
                        updated.write(line)
                updated.write("\n\n## {}\n\n".format(REFERENCE_TITLE))
                for ref in self.references:
                    line = format_ref_to_md(ref)
                    updated.write(line + "\n\n")
                    if self.verbose:
                        print(line)
        os.rename(tmp_file, self.original_file)
                    
            
                    

# !SECTION! Main program

def main_cli():
    t = TestTangler(argv[1], verbose=False)
    t.process()

    
if __name__ == "__main__":  main_cli()
    
