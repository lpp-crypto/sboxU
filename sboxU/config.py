
# !SECTION! Threads number

"""The maximum number of threads used in each function. """
MAX_N_THREADS = int(8)

def n_threads_from_sbox_size(n):
    """Customize this function to adjust the number of threads used by multi-threaded functions.

    The idea is to prevent the use of many threads when each of them actually has very little to do.

    Args:
        - n: the bitlength of an S-box input.

    Returns:
        The number of threads to use. The higher the input, the higher the number of threads.
    """
    if n <= 6:
        return 1
    elif n <= 8:
        return 2
    else:
        return MAX_N_THREADS



# !SECTION! Constants

"""The number of bits of precision to use when performing high precision real arithmetic.

"""
DEFAULT_HIGH_PRECISION = 120


# !SECTION! Terminal output appearance 

SECTION_TEMPLATES = [
    "\n[bold underline]{}[/bold underline]\n\n",                   # title
    "\n[bold underline purple]{:2d} {}[/bold underline purple]\n", # section
    "[blue underline]{:2d}.{:2d}) {} [/blue underline]"            # subsection
]

KEYWORD_TEMPLATES = {
    "fail" : "[bold red]{}[/bold red]\n",
    "success" : "[bold green]{}[/bold green]\n",
}
