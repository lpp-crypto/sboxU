from sboxUv2.core import \
    Sb, \
    degree_spectrum, algebraic_degree

from sboxUv2.statistics import \
    differential_spectrum, \
    walsh_spectrum, absolute_walsh_spectrum

from sboxUv2.ccz import \
    thickness_spectrum, \
    get_WalshZeroesSpaces, \
    ccz_equivalent_function, extended_affine_equivalences

from sboxUv2.apn import \
    get_WalshZeroesSpaces_quadratic_apn, \
    sigma_multiplicities, \
    apn_ea_mugshot, apn_ea_mugshot_from_spectra

from sboxUv2.databases import *

                 

class APNFunctions(FunctionsDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "apn_functions.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains additional
    logic to handle the specifics of APN functions, and in particular
    of their CCZ-equivalence class. Here, "functions" should be
    thought of much more as "extended affine equivalence class
    representative" rather than function.
    
    This class provides another table containing the bases of all the
    spaces of dimension n contained with the Walsh zeroes of a
    function. In order to both save space and store the structure of a
    CCZ-equivalence, each APN function is stored along with the
    identifier of the Walsh spaces of its CCZ-equivalence class, and
    the FastLinearMapping that must be applied to it to obtain its own
    Walsh Zeroes.

    """
    
    def __init__(self, db_file):
        # !IDEA! have a max_degree and a min_degree?
        super().__init__(
            db_file,
            {
                "lut" : "BLOB", # the lookup table of the representative
                "n" : "INTEGER",
                "m" : "INTEGER",
                "linearity" : "INTEGER",
                "thickness" : "INTEGER",
                "degree" : "INTEGER",
                "bibliography" : "TEXT",
                "ccz_id": "INTEGER",
                "mugshot" : "BLOB"
            }
        )
        try:
            self.cursor.execute("SELECT COUNT(ccz_id) FROM {}".format(self.functions_table))
            self.number_of_ccz_classes = self.cursor.fetchall()[0][0]
        except:
            self.number_of_ccz_classes = 0
                
    
    def __str__(self):
        return "APN function DB containing {} EA-classes from {} CCZ-classes".format(
            self.number_of_functions,
            self.number_of_ccz_classes
        )



    # def insert_new_ea_repr(self, s, bibliography):
    #     sb = Sb(s)
    #     differential_spec = differential_spectrum(sb)
    #     if differential_spec.maximum() != 2:
    #         raise Exception("Trying to add a non-APN function to the APN function database: {}".format(lut))
    #     # spectra
    #     walsh_spec = walsh_spectrum(sb)
    #     lin = walsh_spec.absolute().maximum()
    #     degree_spec = degree_spectrum(sb)
    #     deg = degree_spec.maximum()
    #     thk_spec = thickness_spectrum(sb)
    #     thk = thk_spec.maximum()
    #     sig_spec = sigma_multiplicities(sb, k=4)
    #     # we assume that it is from a 
    #     self.number_of_ccz_classes += 1
    #     if deg == 2:
    #         mug = apn_ea_mugshot(sb)
    #     else:
    #         mug = apn_ea_mugshot_from_spectra(walsh_spec,
    #                                           degree_spec,
    #                                           sig_spec,
    #                                           thk_spec)
    #     # inserting the function 
    #     to_insert = {
    #         "lut" : sb.to_bytes(),
    #         "n" : sb.get_input_length(),
    #         "m" : sb.get_output_length(),
    #         "bibliography" : bibliography,
    #         "linearity" : lin,
    #         "degree" : deg,
    #         "thickness" : thk,
    #         "ccz_id" : self.number_of_ccz_classes,
    #         "mugshot" : mug
    #     }
    #     return self.insert_function(to_insert)


    def insert_full_ccz_equivalence_class(self, s, bibliography):
        sb = Sb(s)
        differential_spec = differential_spectrum(sb)
        if differential_spec.maximum() != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(sb))
        encoded = sb.to_bytes()
        # linear
        abs_walsh_spec = absolute_walsh_spectrum(sb)
        lin = abs_walsh_spec.maximum()
        inserted_ids = []
        if algebraic_degree(sb) == 2: # if the function is quadratic,
                                       # we compute automorphisms
                                       # inserting the spaces
            ws = get_WalshZeroesSpaces_quadratic_apn(sb)
            quadratic = True
        else:
            ws = get_WalshZeroesSpaces(sb)
            quadratic = False
        # inserting all the functions
        for L in ws.get_mappings():
            new_sb = ccz_equivalent_function(sb, L)
            new_ws = ws.image_by(L.transpose().inverse())
            new_thk_spec = new_ws.thickness_spectrum()
            new_degree_spec = degree_spectrum(new_sb)
            new_sigma_mult = sigma_multiplicities(new_sb, k=4)
            if quadratic:
                worth_adding = True
            else:
                worth_adding = self.is_new(
                    new_sb,
                    degree_spec=new_degree_spec,
                    abs_walsh_spec=abs_walsh_spec,
                    sigma_mult=new_sigma_mult,
                    thk_spec=new_thk_spec,
                    ccz_id=self.number_of_ccz_classes
                )
            if worth_adding:
                if new_degree_spec.maximum() == 2:
                    mug = apn_ea_mugshot(new_sb)
                else:
                    mug = apn_ea_mugshot_from_spectra(
                        abs_walsh_spec,
                        new_degree_spec,
                        new_sigma_mult,
                        new_thk_spec
                    )
                to_insert = {
                    "lut" : new_sb.to_bytes(),
                    "n" : new_sb.get_input_length(),
                    "m" : new_sb.get_output_length(),
                    "bibliography" : bibliography,
                    "linearity" : lin,
                    "degree" : new_degree_spec.maximum(),
                    "thickness" : new_thk_spec.maximum(),
                    "ccz_id" : self.number_of_ccz_classes,
                    "mugshot" : mug
                }
                inserted_ids.append(self.insert_function(to_insert))
        self.number_of_ccz_classes += 1
        return inserted_ids
    

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        entry["sbox"] = Sb(entry["lut"])
        return entry

    
    def is_new(self,
                   s,
                   degree_spec=None,
                   abs_walsh_spec=None,
                   sigma_mult=None,
                   thk_spec=None,
                   ccz_id=None
        ):
        """
        # !TODO! docstring 
        """
        sb = Sb(s)
        if degree_spec == None:
            degree_spec = degree_spectrum(sb)
        if degree_spec.maximum() == 2:
            mug = apn_ea_mugshot(sb)
        else:
            if abs_walsh_spec == None:
                abs_walsh_spec = absolute_walsh_spectrum(sb)
            if sigma_mult == None:
                sigma_mult = sigma_multiplicities(sb)
            if thk_spec == None:
                thk_spec = thickness_spectrum(sb)
            mug = apn_ea_mugshot_from_spectra(
                abs_walsh_spec,
                degree_spec,
                sigma_mult,
                thk_spec
            )
        query = {"mugshot" : mug}
        # query = {}
        if ccz_id != None:
            query["ccz_id"] = ccz_id
        candidates = self.query_functions(query)
        if candidates == []:
            return True
        else:
            # checking all the functions with a similar mugshot
            for entry in candidates:
                if sb == entry["sbox"]:
                    return False
                test_result = extended_affine_equivalences(
                    sb.lut(),
                    entry["sbox"].lut(),
                    single_non_trivial_answer=False,
                    )
                if len(test_result) > 0:
                    return False
            return True
