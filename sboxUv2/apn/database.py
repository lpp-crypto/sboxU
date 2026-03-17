from sboxUv2.core import \
    Sb, \
    degree_spectrum, algebraic_degree, quadratic_compact_representation, quadratic_sbox_from_compact_representation

from sboxUv2.statistics import \
    differential_spectrum, \
    walsh_spectrum, absolute_walsh_spectrum, linearity

from sboxUv2.ccz import \
    thickness_spectrum, \
    get_WalshZeroesSpaces, \
    ccz_equivalent_function, \
    are_ea_equivalent, \
    are_ccz_equivalent\
   

from sboxUv2.apn import \
    get_WalshZeroesSpaces_quadratic_apn, \
    sigma_multiplicities, \
    apn_ea_mugshot, apn_ea_mugshot_from_spectra, ccz_equivalent_quadratic_function, ea_mappings_from_ortho_derivative
    

from sboxUv2.databases import *
import hashlib

                 

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
        if self.new_db:
            self.number_of_ccz_classes = 0
        else:
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
            raise Exception("Trying to add a non-APN function to the APN function database: \nspec={}\ns={}".format(differential_spec, sb))
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
            new_L = L.transpose()
            new_L = new_L.inverse()
            new_ws = ws.image_by(new_L)
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
        if ccz_id != None:
            query["ccz_id"] = ccz_id
        
        candidates = self.query_functions(query)
        if candidates == []:
            print("+ no candidate")
            return True
        else:
            # checking all the functions with a similar mugshot
            print("\n\n", len(candidates))
            print(mug)
            print(sb)
            print("___")
            for entry in candidates:
                print(entry["id"], entry["ccz_id"])
                print(entry["mugshot"])
                if sb == entry["sbox"]:
                    print("Identic to {}".format(entry["id"]))
                    return False
                
                if are_ea_equivalent(sb.lut(), entry["sbox"].lut()):
                    print("EA to entry {}".format(entry["id"]))
                    print(entry["sbox"])
                    return False
            return True


class APNQuadraticFunctions_ccz_only(FunctionsDB):
    """This class is expected to be bundled with a literal TinySQL
    database file called "apn_functions.db", and allows an easy
    interaction with it.

    It builds upon the `FunctionDB` class, and contains additional
    logic to handle the specifics of APN functions, and in particular
    of their CCZ-equivalence class. Here, "functions" should be
    thought of much more as "extended affine equivalence class
    representative" rather than function.
    
    This class is a compact version of APNFunctions to use when there are
    too many functions when we include EA-classes. 

    """

  
    def __init__(self, db_file):
        # !IDEA! have a max_degree and a min_degree?
        super().__init__(
            db_file,
            {
                "qcr" : "BLOB",
                "mugshot" : "BLOB",
                "linearity" : "INTEGER"
            }
        )
        if self.new_db:
            self.number_of_ccz_classes = 0
        else:
            try:
                self.cursor.execute("SELECT COUNT(id) FROM {}".format(self.functions_table))
                self.number_of_ccz_classes = self.cursor.fetchall()[0][0]
            except:
                self.number_of_ccz_classes = 0
                
    
    def __str__(self):
        return "APN function DB containing {} CCZ-classes".format(
            self.number_of_functions

        )


    def insert_quadratic_ccz_representative(self, s):
        """
        In this version, we only insert a ccz representative
        """
        
        sb = Sb(s)
        differential_spec = differential_spectrum(sb)
        if differential_spec.maximum() != 2:
            raise Exception("Trying to add a non-APN function to the APN function database: {}".format(sb))

        if algebraic_degree(sb) == 2: 
            mug = apn_ea_mugshot(sb)
            # We hash the mugshot for memory concerns
            h = hashlib.sha256()
            h.update(mug)
            mug = h.digest()

            worth_adding = self.is_new(sb,mug)
        else:
            s_quad = ccz_equivalent_quadratic_function(sb.lut())
            if s_quad != []:
                sb = s_quad
                mug = apn_ea_mugshot(sb)
                # We hash the mugshot for memory concerns
                h = hashlib.sha256()
                h.update(mug)
                mug = h.digest()
                worth_adding = self.is_new(sb,mug)
            else:
                raise Exception("Trying to add a non ccz_quadratic APN function to the database: {}".format(sb))



        if worth_adding:
            to_insert = {
                "qcr" : bytearray(quadratic_compact_representation(sb.lut())),
                "mugshot" : mug,
                "linearity" : linearity(sb.lut())
            }
            inserted_id = self.insert_function(to_insert)
            self.number_of_ccz_classes += 1
            return inserted_id
        return None
    

        
    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        # post-processing
        #entry["sbox"] = Sb(quadratic_sbox_from_compact_representation(entry["qcr"],8,8))
        return entry

    
    def is_new(self,
                   s,
                   mug=None
        ):
        """
        # !TODO! docstring 
        """
        sb = Sb(s)
        
        if mug == None:
            # Computing the mugshot of s depending on its degree
            if degree_spec == None:
                degree_spec = degree_spectrum(sb)
            if degree_spec.maximum() == 2:
                mug = apn_ea_mugshot(sb)
            else:
                # If s is ccz quadratic, we work with a quadratic representative
                s_quad = ccz_equivalent_quadratic_function(sb)
                if s_quad == []:
                    # !! TODO !!
                    # Decide if we put more
                    mug = absolute_walsh_spectrum(sb)
                else:
                    mug = apn_ea_mugshot(Sb(s_quad))


            # We hash the mugshot for memory concerns
            h = hashlib.sha256()
            h.update(mug)
            mug = h.digest()


        query = {"mugshot" : mug}
        candidates = self.query_functions(query)
        if candidates == []:
            #print("- New Mugshot")
            return True
        else:
            for entry in candidates:
                print("Same Mugshot as Function {}".format(entry["id"]))
                print(entry["mugshot"])
                
                # Trivial Case
                if quadratic_compact_representation(sb.lut()) == entry["qcr"]:
                    return False
                
                # Testing EA-equivalence
                qcr = entry["qcr"]
                lut = quadratic_sbox_from_compact_representation(qcr,sb.get_input_length(),sb.get_output_length())
                mappings = ea_mappings_from_ortho_derivative(sb.lut(), lut)
                if mappings != []:
                    print("Function EA equivalent to  id = {} in the database".format(entry["id"]))
                    return False

            return True

    def insert_many_functions(self,entries):
        
        start_id = self.number_of_functions
        end_id  = start_id
        to_insert = []
        for e in entries:
            e["id"] = end_id
            end_id +=1
            e_tuple = tuple([e[k] for k in sorted(self.row_structure.keys())])
            to_insert.append(e_tuple)
        try:
            #self.cursor.executemany(self.function_insertion_query,to_insert)
            self.cursor.executemany( "INSERT INTO functions(id, linearity, mugshot, qcr) VALUES (?,?,?,?)",to_insert)
            self.number_of_functions = end_id
            return list(range(start_id,end_id+1))
        except:
            raise Exception("Insertion of many entries failed \n")

