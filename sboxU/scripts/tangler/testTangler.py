from sys import argv


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
        
        
    def finalize_current_block(self) -> None:            
        local_indent = 2*self.INDENT if self.experiment_has_started else ""
        self.tangled.write(local_indent + self.BLOCK_BEGIN)
        for line in self.current_block:
            self.tangled.write(local_indent + line)
        self.tangled.write(local_indent + self.BLOCK_END)
        self.current_block = None    
        if self.verbose:
            print("tangled {} lines".format(len(self.current_block)))

            
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
        
        
    def absorb_line(self, line : str):
        if line[0:3] == "```":
            if self.current_block != None:
                self.finalize_current_block()
            else:
                self.current_block = []
        elif isinstance(self.current_block, (list)):
            self.current_block.append(line)
        elif line[0] == "#":
            self.process_section(line)
        

    def process(self):
        if self.verbose:
            print("Tangling from {}".format(self.original_file))
        with open(self.original_file, "r") as original:
            self.original = original
            with open(self.tangled_file, "w") as tangled:
                self.tangled = tangled
                self.tangled.write(self.HEADER)
                for line in original.readlines():
                    if len(line) > 3:
                        self.absorb_line(line)
                self.tangled.write(self.INDENT + "return exit_code()\n\n")
                self.tangled.write(self.FOOTER)
                    

# !SECTION! Main program

def main_cli():
    t = TestTangler(argv[1])
    t.process()

    
if __name__ == "__main__":  main_cli()
    
