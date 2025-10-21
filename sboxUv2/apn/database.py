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



    def insert_new_ea_repr(self, s, bibliography):
        sb = Sb(s)
        differential_spec = differential_spectrum(sb)
        if differential_spec.maximum() != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(lut))
        encoded = sb.to_bytes()
        # spectra
        walsh_spec = walsh_spectrum(sb)
        lin = walsh_spec.absolute().maximum()
        degree_spec = degree_spectrum(sb)
        deg = degree_spec.maximum()
        thk_spec = thickness_spectrum(sb)
        thk = thk_spec.maximum()
        sig_spec = sigma_multiplicities(sb)
        # we assume that it is from a 
        self.number_of_ccz_classes += 1
        # inserting the function 
        to_insert = {
            "lut" : encoded,
            "n" : n,
            "m" : m,
            "bibliography" : bibliography,
            "linearity" : lin,
            "degree" : deg,
            "thickness" : thk,
            "ccz_id" : self.number_of_ccz_classes,
            "mugshot" : apn_ea_mugshot_from_spectra(walsh_spec,
                                                    degree_spec,
                                                    sig_spec,
                                                    thk_spec)
        }
        return self.insert_function(to_insert)


    def insert_full_ccz_equivalence_class(self, s, bibliography):
        sb = Sb(s)
        differential_spec = differential_spectrum(lut)
        if differential_spec.maximum() != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(lut))
        encoded = sb.to_bytes()
        # linear
        walsh_spec = walsh_spectrum(sb)
        lin = walsh_spec.absolute().maximum()
        # !CONTINUE! stopped here. 
        if algebraic_degree(lut) == 2: # if the function is quadratic,
                                       # we compute automorphisms
                                       # inserting the spaces
            quadratic_case = True
            autom = [FastLinearMapping(L)
                     for L in graph_automorphisms_of_apn_quadratic(lut)]
            Ls = spaces.reduced_linear_mappings(autom)
        else:
            quadratic_case = False
            Ls = spaces.Ls
        self.cursor.execute(self.spaces_insertion_query, (
            self.number_of_ccz_classes,
            spaces.to_blob()
        ))
        # inserting all the functions
        for L in Ls:
            new_lut = apply_mapping_to_graph(lut, L)
            if quadratic:
                worth_adding = True
            else:
                worth_adding = (self.is_present(new_lut, test_ccz=True)[0] == "absent")
            if worth_adding:
                L_inv_T = L.transpose().inverse()
                new_spaces = L_inv_T * spaces
                new_thk_spec = new_spaces.thickness_spectrum()
                t1 = pretty_spectrum(thickness_spectrum(new_lut))
                t2 = pretty_spectrum(new_thk_spec)
                if t1 != t2:
                    FAIL("well there's your problem: {} != {}".format(t1, t2))
                new_degree_spec =degree_spectrum(new_lut)
                to_insert = {
                    "lut" : encode_lut(new_lut, n),
                    "n" : n,
                    "m" : m,
                    "bibliography" : bibliography,
                    "linearity" : lin,
                    "degree" : max(new_degree_spec.keys()),
                    "thickness" : max(new_thk_spec.keys()),
                    "ccz_id" : self.number_of_ccz_classes,
                    "walsh_L" : encode_lut(L_inv_T.masks, n+m),
                    "mugshot" : mugshot(new_lut,
                                        walsh_spec=walsh_spec,
                                        degree_spec=new_degree_spec,
                                        differential_spec=differential_spec,
                                        thickness_spec=new_thk_spec)
                }
                self.insert_function(to_insert)
        self.number_of_ccz_classes += 1
        return self.number_of_ccz_classes
    

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        entry["lut"] = decode_lut(entry["lut"])
        W = self.fetch_WalshZeroesSpaces(entry["ccz_id"])
        L = FastLinearMapping(decode_lut(entry["walsh_L"]))
        entry["walsh_spaces"] = L * W
        return entry

    
    def is_present(self,
                   s,
                   degree_spec=None,
                   walsh_spec=None,
                   differential_spec=None,
                   thickness_spec=None,
                   test_ccz=True):
        """Returns one of three things:
        
        - `["present", index]` if there is already a function
          extended-affine equivalent to `s` in the data-base where
          `index` is its `id`,
        
        - `["absent"]` if there isn't, and
        
        - [`"maybe"`, index]` if at least one function with a similar
          mugshot is present, that function having an `id` equal to
          `index`.

        If `test_ccz` is set to `True` (the default case), then a
        CCZ-equivalence test is performed using SAGE's implementation
        of code equivalence to get a definitive answer. `"maybe"` can
        only be returned if `test_ccz` is set to `False`.

        """
        n, m = get_block_lengths(s)
        encoded = encode_lut(s, n)
        mug = mugshot(s,
                      degree_spec=degree_spec,
                      walsh_spec=walsh_spec,
                      differential_spec=differential_spec,
                      thickness_spec=thickness_spec)
        candidates = self.query_functions({"mugshot" : mug})
        if candidates == []:
            return ["absent", 0]
        elif test_ccz:
            # checking all the functions with a similar mugshot
            for entry in candidates:
                if entry["lut"] == encoded:
                    return ["present", entry["id"], entry["ccz_id"]]
                else:
                    #print("manually testing CCZ equivalence")
                    if are_ccz_equivalent(entry["lut"], s):
                        return ["present", entry["id"], entry["ccz_id"]]
            # if we checked for CCZ-equivalence and didn't find any
            # CCZ-equivalent function, then `s` is new...
            return ["absent", 1]
        else:
            # ... but if we didn't check, then we can't be sure.
            return ["maybe", candidates[0]["id"]]
