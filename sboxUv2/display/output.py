
# !SECTION! General setup


# !SUBSECTION! Importing classes needed for tests

from collections import defaultdict

# !SUBSECTION! Importing pretty printing utilities from Rich 

from rich.console import Console
from rich.theme import Theme


# !SUBSECTION! Importing tools to run a timer 

import datetime
import time

from math import floor


# !SUBSECTION! Setting up global variables

from ..config import SECTION_TEMPLATES
CONSOLE = Console(theme=Theme({}, inherit=False))
ONGOING_EXPERIMENT = None




# !SECTION! The main pretty printing function 

def pprint(*args):
    result = ""
    all_pretty = True
    for x in args:
        if hasattr(x, "__rich_str__"):
            result += x.__rich_str__()
        elif isinstance(x, (dict, defaultdict)):
            result = "{ "
            for k in sorted(x.keys()):
                result += "[bold]{}[/bold]: {}, ".format(k, x[k])
            result = result[:-2] + " }"
        else:
            result += str(x)
            all_pretty = False
        result += ", "
    result = result[:-2] # ditching superfluous ", " 
    if all_pretty:
        CONSOLE.print(result)
    else:
        print(result)


        
# !SECTION! Handling sections in the output


# !SUBSECTION! The Timer class


class Chronograph:
    def __init__(self, title):
        self.title = title
        self.start_time = datetime.datetime.now()

    def __str__(self):
        elapsed_time = datetime.datetime.now() - self.start_time
        tot_secs = floor(elapsed_time.total_seconds())
        days = floor(tot_secs / 86400)
        hours = floor((tot_secs % 86400) / 3600)
        minutes = floor((tot_secs % 3600) / 60)
        seconds = (tot_secs % 60) + elapsed_time.total_seconds() - tot_secs
        return "\"{}\" lasted {}s ({})".format(
            self.title,
            elapsed_time.total_seconds(),
            "{:d}d {:02d}h {:02d}m {:5.03f}s".format(
                days,
                hours,
                minutes,
                seconds
        ))


    def __rich_str__(self):
        elapsed_time = datetime.datetime.now() - self.start_time
        tot_secs = floor(elapsed_time.total_seconds())
        days = floor(tot_secs / 86400)
        hours = floor((tot_secs % 86400) / 3600)
        minutes = floor((tot_secs % 3600) / 60)
        seconds = (tot_secs % 60) + elapsed_time.total_seconds() - tot_secs
        return "[blue]{}[/blue] lasted [bold]{}[/bold]s [gray]({})[/gray]".format(
            self.title,
            elapsed_time.total_seconds(),
            "{:d}d {:02d}h {:02d}m {:5.03f}s".format(
                days,
                hours,
                minutes,
                seconds
        ))



# !SUBSECTION!  Starting and finishing sections

class Experiment:
    """The purpose of this class is to simplify the task of generating a pretty output in a terminal that contains more automatically derived information.

    To use it, simply rely on the `with` statement:

    Example:

        with Experiment("your experiment title"):
            section("the first part")
            <your code>
            section("the second part (which is bigger)")
            subsection("the first part of the second part")
            <your code>

    When running your code, the title of the experiment will be displayed, and the sections titles will be highlighted for ease of reading.

    The runtime of each section and subsection is also displayed at the end of each of these.
    
    """
    def __init__(self, title):
        self.title = title
        self.sections_counters = [0]
        self.section_timer = None
        CONSOLE.print(SECTION_TEMPLATES[0].format(self.title))
        self.global_timer = Chronograph("[bold]{}[/bold]".format(self.title))

        
    def section(self, title):
        if self.sections_counters[0] > 0:
            pprint(self.section_timer)
        self.sections_counters = [self.sections_counters[0] + 1]
        self.section_timer = Chronograph(title)
        CONSOLE.print(SECTION_TEMPLATES[1].format(
            self.sections_counters[0],
            title))

        
    def subsection(self, title):
        if len(self.sections_counters) == 1:
            self.sections_counters.append(1)
        else:
            self.sections_counters = [
                self.sections_counters[0],
                self.sections_counters[1] + 1
            ]
        CONSOLE.print(SECTION_TEMPLATES[2].format(
            self.sections_counters[0],
            self.sections_counters[1],
            title
        ))

    def __enter__(self):
        global ONGOING_EXPERIMENT
        ONGOING_EXPERIMENT = self

    
    def __exit__(self, exception_type, exception_value, traceback):
        if self.sections_counters[0] > 0:
            pprint(self.section_timer)
        pprint(self.global_timer)
        

def section(title):
    global ONGOING_EXPERIMENT
    if ONGOING_EXPERIMENT == None:
        raise Exception("sections can only be used within an Experiment")
    ONGOING_EXPERIMENT.section(title)

    
def subsection(title):
    global ONGOING_EXPERIMENT
    if ONGOING_EXPERIMENT == None:
        raise Exception("subsections can only be used within an Experiment")
    ONGOING_EXPERIMENT.subsection(title)
