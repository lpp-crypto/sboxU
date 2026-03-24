from sboxU import *

if __name__ == "__main__":
    print(cite("AC:BonPerTia19"))
    for entry in who_to_cite(extract_bases):
        print("---")
        print(cite_as(entry, "rst"))
    # print("_______________\n\n")
    # print("")
    # # print_full_biblio_rst()
