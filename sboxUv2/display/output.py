from rich.console import Console
from rich.theme import Theme

console = Console(theme=Theme({}, inherit=False))



def pprint(*args):
    result = ""
    all_pretty = True
    for x in args:
        if hasattr(x, "__rich_str__"):
            result += x.__rich_str__()
        else:
            result += str(x)
            all_pretty = False
        result += ", "
    result = result[:-2] # ditching superfluous ", " 
    if all_pretty:
        console.print(result)
    else:
        print(result)
