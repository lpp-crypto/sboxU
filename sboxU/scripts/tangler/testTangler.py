from sys import argv
import re
import os

from sboxU.biblio import *

REFERENCE_HEADER = "## References"
BLOCK_DELIMITER = "```"


# !SECTION! The TestTangler class

class TestTangler:
    """This class is used to parse a mardown file containing SAGE code with some minimal structure requirements, and to turn it into an sboxU-style `Experiment`.
    
    """

    PREAMBLE_TITLE  = "Preamble"
    PREAMBLE  = "# Tangled from {}\n\n{}\n\nwith Experiment('{}'):\n"
    INDENT    = "    "
    BLOCK_BEGIN = "# --- { \n"
    BLOCK_END   = "# --- } \n"
    HEADER = "import sys\nfrom sage.all import *\nfrom sboxU import *\n"
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
        
        
    def finalize_current_block(self) -> None:            
        local_indent = 2*self.INDENT if self.experiment_has_started else ""
        self.tangled.write(local_indent + self.BLOCK_BEGIN)
        for line in self.current_block:
            self.tangled.write(local_indent + line)
        self.tangled.write(local_indent + self.BLOCK_END)
        if self.verbose:
            print("tangled {} lines".format(len(self.current_block)))
        self.current_block = None    

            
    def start_experiment(self, title : str) -> None:
        self.tangled.write("def main_test():\n")
        self.tangled.write(self.INDENT + "with Experiment('{}'):\n".format(title))
        self.experiment_has_started = True
        if self.verbose:
            print("Experiment: ", title)
        
            
    def process_section(self, line : str) -> None:
        title = line.split("#")[-1][:-1] # dropping final new line
        if line[:2] == "# ":
            self.experiment_title = title
        elif line[:3] == "## ":
            if self.PREAMBLE_TITLE not in title:
                if self.experiment_has_started == False:
                    self.start_experiment(title)    
                self.tangled.write(2*self.INDENT + "section('{}')\n".format(title))
        elif line[:4] == "### ":
            if self.PREAMBLE not in title:
                self.tangled.write(2*self.INDENT + "subsection('{}')\n".format(title))
        elif line[:4] == "####": # all other sections
            if self.PREAMBLE not in title:
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
                self.tangled.write(self.HEADER)
                for line in original.readlines():
                    if len(line) > 1:
                        self.absorb_line(line)
                        if line.startswith(REFERENCE_HEADER):
                            break
                self.tangled.write(self.INDENT + "return exit_code()\n\n")
                self.tangled.write(self.FOOTER)
        if len(self.references) > 0:
            print("handling refs")
            self.add_references()

            
    def add_references(self) -> None:
        tmp_file = self.original_file + ".tmp"
        with open(self.original_file, "r") as original:
            with open(tmp_file, "w") as updated:
                for line in original.readlines():
                    if line.startswith(REFERENCE_HEADER):
                        break
                    else:
                        updated.write(line)
                updated.write("\n\n{}\n\n".format(REFERENCE_HEADER))
                for ref in self.references:
                    line = format_ref_to_md(ref)
                    updated.write(line + "\n\n")
                    print(line)
        os.rename(tmp_file, self.original_file)
                    
            
                    

# !SECTION! Main program

def main_cli():
    t = TestTangler(argv[1], verbose=False)
    t.process()

    
if __name__ == "__main__":  main_cli()
    
