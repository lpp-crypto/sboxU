def stylize(line, desc):
    if "b" in desc:
        line = "\033[1m" + line + "\033[0m"
    if "r" in desc:
        line = "\033[91m" + line + "\033[0m"
    return line
