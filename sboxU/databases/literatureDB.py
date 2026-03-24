from sboxU.core import S_box, get_sbox, degree_spectrum
from sboxU.statistics import differential_spectrum, absolute_walsh_spectrum


from sboxU.databases.database import FunctionsDB


class LiteratureSBoxes(FunctionsDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "literature_sboxes.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains just a bit of
    logic on top to handle the computation of all the S-box properties
    we are interested in.

    # !TODO! add representative of linear equivalence classes; and the
    # !logic to use it

    """
    
    def __init__(self):
        super().__init__(
            "literature_sboxes.db",
            {
                "lut" : "BLOB",
                "n" : "INTEGER",
                "m" : "INTEGER",
                "differential_uniformity" : "INTEGER",
                "linearity" : "INTEGER",
                "degree" : "INTEGER",
                "bibliography" : "TEXT",
                "name" : "TEXT",
                "usage": "TEXT"
            }
        )

        
    def insert_function_from_lut(self, s, name, bibliography, usage=None):
        sb = get_sbox(s)
        n, m = sb.get_input_length(), sb.get_output_length()
        encoded = sb.to_bytes()
        # differential
        diff_spec = differential_spectrum(sb)
        diff_unif = diff_spec.maximum()
        # linear
        walsh_spec = absolute_walsh_spectrum(sb)
        lin = walsh_spec.maximum()
        # degree
        deg = degree_spectrum(sb).maximum()
        # usage
        if usage == None:
            usage = name
        # building the query
        to_insert = {
            "lut" : encoded,
            "n" : n,
            "m" : m,
            "bibliography" : bibliography,
            "differential_uniformity" : diff_unif,
            "linearity" : lin,
            "degree" : deg,
            "name" : name,
            "usage" : usage
        }
        self.insert_function(to_insert)

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        entry["sbox"] = get_sbox(entry["lut"], name=entry["name"])
        return entry
