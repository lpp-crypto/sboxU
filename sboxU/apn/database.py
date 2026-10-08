from os.path import dirname

from sboxU.core import \
    get_sbox, get_F2AffineMap, \
    degree_spectrum, algebraic_degree, quadratic_compact_representation, quadratic_sbox_from_compact_representation

from sboxU.statistics import \
    differential_spectrum, \
    walsh_spectrum, absolute_walsh_spectrum, linearity

from sboxU.ccz import \
    thickness_spectrum, \
    get_WalshZeroesSpaces, \
    ccz_equivalent_function, \
    ccz_equivalences, \
    are_ea_equivalent, \
    are_ea_equivalent_from_vq, \
    are_ccz_equivalent

from sboxU.apn import \
    get_WalshZeroesSpaces_quadratic_apn, \
    sigma_multiplicities, \
    apn_ea_mugshot, apn_ea_mugshot_from_spectra, ccz_equivalent_quadratic_function, ea_mappings_from_ortho_derivative, \
    automorphisms_from_ortho_derivative


from sboxU.databases import *
import hashlib


def sixBitAPNs():
    """Returns the path to the database of 6-bit functions bundled with sboxU.

    To be used when building an `APNFunctions` object:

    `with APNFunctions(sixBitAPNs()) as db:`
    """
    return dirname(__file__) + "/../scripts/apnDB/apn6.db"



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
            except Exception:
                self.number_of_ccz_classes = 0


    def __str__(self):
        return "APN function DB containing {} EA-classes from {} CCZ-classes".format(
            self.number_of_functions,
            self.number_of_ccz_classes
        )



    # def insert_new_ea_repr(self, s, bibliography):
    #     sb = get_sbox(s)
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
        sb = get_sbox(s)
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
        entry["sbox"] = get_sbox(entry["lut"])
        return entry


    def is_new(self,
                   s,
                   degree_spec=None,
                   abs_walsh_spec=None,
                   sigma_mult=None,
                   thk_spec=None,
                   ccz_id=None
        ):
        """Returns True if and only if `s` is not EA-equivalent to any function already in the database.

        Computes an EA-invariant "mugshot" for `s` and queries the database for matching entries.
        If any candidate is found, EA-equivalence is verified explicitly. Pre-computed spectra can
        be supplied to avoid redundant computation.

        Args:
            s: An S-boxable APN function.
            degree_spec: Pre-computed degree spectrum of `s`, or None to compute it.
            abs_walsh_spec: Pre-computed absolute Walsh spectrum of `s`, or None to compute it.
            sigma_mult: Pre-computed sigma multiplicities of `s`, or None to compute them.
            thk_spec: Pre-computed thickness spectrum of `s`, or None to compute it.
            ccz_id: If given, restrict the search to functions with this CCZ class identifier.

        Returns:
            bool: True if `s` is new (not EA-equivalent to any stored function), False otherwise.
        """
        sb = get_sbox(s)
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
            except Exception:
                self.number_of_ccz_classes = 0


    def __str__(self):
        return "APN function DB containing {} CCZ-classes".format(
            self.number_of_functions

        )


    def insert_quadratic_ccz_representative(self, s):
        """
        In this version, we only insert a ccz representative
        """

        sb = get_sbox(s)
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
        #entry["sbox"] = get_sbox(quadratic_sbox_from_compact_representation(entry["qcr"],8,8))
        return entry


    def is_new(self,
                   s,
                   mug=None,
                   degree_spec=None
        ):
        """Returns True if and only if `s` is not EA-equivalent to any function already in the database.

        This compact version uses a hashed mugshot for fast lookup, relying on the quadratic compact
        representation and EA-mappings from the ortho-derivative for confirmation.

        Args:
            s: An S-boxable APN function.
            mug: A pre-computed (and hashed) mugshot bytes object, or None to compute it.
            degree_spec: Pre-computed degree spectrum of `s`, or None to compute it.

        Returns:
            bool: True if `s` is new (not EA-equivalent to any stored function), False otherwise.
        """
        sb = get_sbox(s)

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
                    mug = apn_ea_mugshot(get_sbox(s_quad))


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
        except Exception as e:
            raise Exception("Insertion of many entries failed \n") from e


# For CCZ class computation for n=8
# Experimental, be careful of edge cases


# !SECTION! Compact multi-degree APN database


# Number of bytes used for the truncated SHA-256 invariant.
# 12 bytes = 96 bits
_INVARIANT_BYTES = 12

# Entry type constants
TYPE_QUADRATIC    = 0   # quadratic APN; representation = quadratic_compact_representation
TYPE_CCZ_MAPPING  = 1   # CCZ-equivalent to a stored quadratic; representation = image vectors + constant of the admissible mapping (little-endian, ceil(2n/8) bytes each)
TYPE_NONQUADRATIC = 2   # non-CCZ-quadratic APN; representation = LUT bytes


def mapping_to_bytes(mapping, n):
    """Serializes an F2AffineMap on F_2^(2n) as its 2n image vectors followed by its constant.

    Each value is stored as ceil(2n/8) bytes in little-endian order.
    Total size: (2n + 1) * ceil(2n/8) bytes.
    """
    bpv = APNFunctions_compact._bytes_per_val(n)
    result = bytearray()
    for iv in mapping.get_image_vectors():
        result += int(iv).to_bytes(bpv, 'little')
    result += int(mapping.get_cste()).to_bytes(bpv, 'little')
    return result


def bytes_to_mapping(rep, n):
    """Deserializes an F2AffineMap on F_2^(2n) from bytes produced by mapping_to_bytes."""
    bpv = APNFunctions_compact._bytes_per_val(n)
    n_ivs = 2 * n
    ivs = [int.from_bytes(rep[i * bpv:(i + 1) * bpv], 'little') for i in range(n_ivs)]
    cstte = int.from_bytes(rep[n_ivs * bpv:(n_ivs + 1) * bpv], 'little')
    mapping = get_F2AffineMap(ivs)
    if cstte != 0:
        mapping = mapping + cstte
    return mapping


class APNFunctions_compact(FunctionsDB):
    """Compact database of n-bit APN functions supporting three entry types:

        TYPE_QUADRATIC (0)   — quadratic APN, stored as a quadratic compact representation.
        TYPE_CCZ_MAPPING (1) — CCZ-equivalent to a quadratic entry in the same ccz_id class, noted quad_ref,
                               stored as the image vectors followed by the constant of the
                               admissible mapping L (each value in ceil(2n/8) bytes,
                               little-endian) such that ccz_equivalent_function(quad_ref, L)
                               gives this function's graph.
        TYPE_NONQUADRATIC (2)— non-CCZ-quadratic APN; stored as raw LUT bytes.

    All entries carry:
        ccz_id    : integer identifying the CCZ-equivalence class (incremented as classes are added).
        invariant : first 12 bytes (_INVARIANT_BYTES = 96 bits) of the SHA-256 of the EA mugshot,
                    used as a fast filter to check for new functions.
        degree    : algebraic degree of the function.
        linearity : linearity of the function (maximum of the absolute Walsh transform).
        thickness : maximum thickness over all Walsh zero spaces of the function.
        source    : integer index into biblio_dict.from_which_paper(n, source) for provenance.
        ccz_size  : number of EA classes in this function's CCZ class, or -1 if unknown.
                    Only filled in for entries whose CCZ class was built via
                    add_quadratic_ccz_class or populate_quadratic_ccz_class.
        aut_size  : size of Aut(q), the graph automorphism group of the class's quadratic
                    reference, or -1 if unknown (same availability as ccz_size).

    The dimension `n` is a construction-time parameter stored as `self.n`.
    """

    def __init__(self, db_file, n):
        self.n = n
        super().__init__(
            db_file,
            {
                "type"           : "INTEGER",
                "ccz_id"         : "INTEGER",
                "invariant"      : "BLOB",
                "representation" : "BLOB",
                "degree"         : "INTEGER",
                "linearity"      : "INTEGER",
                "thickness"      : "INTEGER",
                "source"         : "INTEGER",
                "ccz_size"       : "INTEGER",
                "aut_size"       : "INTEGER",
            }
        )
        # Speeds up invariant-based duplicate lookups (is_new, same_invariant) from an
        # O(N) full table scan to an O(log N) index lookup.  IF NOT EXISTS makes this
        # a cheap no-op on databases that already have the index, and backfills it on
        # existing databases created before this index was introduced.
        self.cursor.execute(
            "CREATE INDEX IF NOT EXISTS idx_invariant ON {} (invariant)".format(self.functions_table)
        )
        if self.new_db:
            self.number_of_ccz_classes = 0
        else:
            try:
                self.cursor.execute(
                    "SELECT COUNT(DISTINCT ccz_id) FROM {}".format(self.functions_table)
                )
                self.number_of_ccz_classes = self.cursor.fetchone()[0]
            except Exception:
                self.number_of_ccz_classes = 0


    def __str__(self):
        return "APNFunctions_compact (n={}) — {} entries from {} CCZ-classes".format(
            self.n, self.number_of_functions, self.number_of_ccz_classes
        )


    # Since computing the thickness is as hard as computing the Walsh Zeroes
    # We do not compute it by default. The thickness is then defaulted to -1
    # To add the thickness, either use the populate function for quadratic function
    # that will add the entire ccz-class as well, or add it by hand using update database
    _THICKNESS_NOT_COMPUTED = -1


    # !SUBSECTION! Invariant helpers

    @staticmethod
    def _make_invariant(mug):
        """Returns the first _INVARIANT_BYTES bytes of the SHA-256 of `mug`."""
        h = hashlib.sha256()
        h.update(mug)
        return h.digest()[:_INVARIANT_BYTES]


    def _compute_invariant(self, sb):
        """Computes the stored invariant bytes for an S-box, using the appropriate mugshot."""
        if algebraic_degree(sb) == 2:
            mug = apn_ea_mugshot(sb)
        else:
            abs_walsh = absolute_walsh_spectrum(sb)
            deg_spec  = degree_spectrum(sb)
            thk_spec  = thickness_spectrum(sb)
            mug = apn_ea_mugshot_from_spectra(
                abs_walsh, deg_spec, sigma_multiplicities(sb, k=4), thk_spec
            )
        return bytearray(self._make_invariant(mug)), mug


    # !SUBSECTION! Mapping serialization helpers

    @staticmethod
    #!TODO! Check if this function needs to be removed for consistency
    def _bytes_per_val(n):
        """Number of bytes needed to store one value from F_2^(2n)."""
        return (2 * n + 7) // 8

    def _mapping_to_bytes(self, mapping):
        """mapping_to_bytes for this database's n."""
        return mapping_to_bytes(mapping, self.n)

    def _bytes_to_mapping(self, rep):
        """bytes_to_mapping for this database's n."""
        return bytes_to_mapping(rep, self.n)


    # !SUBSECTION! Row parsing

    def parse_function_from_row(self, row):
        entry = {}
        for i, column in enumerate(sorted(self.row_structure.keys())):
            entry[column] = row[i]
        n   = self.n
        t   = entry["type"]
        rep = entry["representation"]
        if t == TYPE_QUADRATIC:
            entry["sbox"] = get_sbox(quadratic_sbox_from_compact_representation(rep, n, n))
        elif t == TYPE_CCZ_MAPPING:
            quad_entries = self.query_functions({"ccz_id": entry["ccz_id"], "type": TYPE_QUADRATIC})
            if not quad_entries:
                raise RuntimeError(
                    "parse_function_from_row: no TYPE_QUADRATIC entry for ccz_id={}".format(
                        entry["ccz_id"]
                    )
                )
            quad_sb = quad_entries[0]["sbox"]
            mapping = self._bytes_to_mapping(rep)
            entry["sbox"] = ccz_equivalent_function(quad_sb, mapping)
        elif t == TYPE_NONQUADRATIC:
            entry["sbox"] = get_sbox(bytes(rep))
        else:
            raise ValueError("parse_function_from_row: unknown entry type {}".format(t))
        return entry


    # !SUBSECTION! Recovering the LUT from an entry

    def get_lut(self, entry):
        """Returns the lookup table of the function described by `entry`."""
        if "sbox" in entry:
            return entry["sbox"].lut()
        return self.parse_function_from_row(
            tuple(entry[c] for c in sorted(self.row_structure.keys()))
        )["sbox"].lut()


    def same_invariant(self, s):
        """Returns the row ids of entries whose invariant matches that of `s`.

        Uses a raw SQL lookup on the invariant column; no function reconstruction is performed.

        Args:
            s: an S-boxable APN function.

        Returns:
            A list of integer row ids (possibly empty).
        """
        sb = get_sbox(s)
        invariant, _ = self._compute_invariant(sb)
        self.cursor.execute(
            "SELECT id FROM {} WHERE invariant = ?".format(self.functions_table),
            (invariant,)
        )
        return [row[0] for row in self.cursor.fetchall()]


    def is_new(self, s, mug=None):
        """Returns True iff `s` is not EA-equivalent to any function already in the database.

        Looks up candidates by their stored invariant.  When `s` is CCZ-quadratic,
        uses the fast Walsh-zero-space EA check (`are_ea_equivalent_from_vq`) against
        candidates that are themselves CCZ-quadratic; otherwise falls back to
        `are_ea_equivalent`.

        Args:
            s:   an S-boxable APN function.
            mug: pre-computed raw EA mugshot (output of `apn_ea_mugshot` /
                 `apn_ea_mugshot_from_spectra`), or None to compute it.

        Returns:
            bool: True if `s` is not EA-equivalent to any stored function.
        """
        sb     = get_sbox(s)
        sb_lut = sb.lut()
        sb_is_ccz_quadratic = len(ccz_equivalent_quadratic_function(sb_lut)) > 0

        if mug is None:
            invariant, mug = self._compute_invariant(sb)
        else:
            invariant = bytearray(self._make_invariant(mug))

        self.cursor.execute(
            "SELECT * FROM {} WHERE invariant = ?".format(self.functions_table),
            (invariant,)
        )
        rows = self.cursor.fetchall()
        if not rows:
            return True

        sorted_columns = sorted(self.row_structure.keys())
        for row in rows:
            entry = {col: row[i] for i, col in enumerate(sorted_columns)}
            candidate_lut = self.get_lut(entry)
            if sb_is_ccz_quadratic and len(ccz_equivalent_quadratic_function(candidate_lut)) > 0:
                if are_ea_equivalent_from_vq(sb_lut, candidate_lut):
                    return False
            else:
                if are_ea_equivalent(sb_lut, candidate_lut):
                    return False

        return True


    # !SUBSECTION! Insertion helpers

    def _insert_entry(self, entry_type, ccz_id, invariant, representation,
                      deg, lin, thk, source, ccz_size=-1, aut_size=-1):
        to_insert = {
            "type"           : int(entry_type),
            "ccz_id"         : int(ccz_id),
            "invariant"      : bytearray(invariant),
            "representation" : bytearray(representation),
            "degree"         : int(deg),
            "linearity"      : int(lin),
            "thickness"      : int(thk),
            "source"         : int(source),
            "ccz_size"       : int(ccz_size),
            "aut_size"       : int(aut_size),
        }
        return self.insert_function(to_insert)


    def _insert_entries(self, rows):
        """Batch analog of `_insert_entry`.

        Args:
            rows: a list of (entry_type, ccz_id, invariant, representation,
                  deg, lin, thk, source, ccz_size, aut_size) tuples.

        Returns:
            A list of row ids, one per element of `rows`, in the same order.
        """
        to_insert = [
            {
                "type"           : int(entry_type),
                "ccz_id"         : int(ccz_id),
                "invariant"      : bytearray(invariant),
                "representation" : bytearray(representation),
                "degree"         : int(deg),
                "linearity"      : int(lin),
                "thickness"      : int(thk),
                "source"         : int(source),
                "ccz_size"       : int(ccz_size),
                "aut_size"       : int(aut_size),
            }
            for entry_type, ccz_id, invariant, representation, deg, lin, thk, source, ccz_size, aut_size in rows
        ]
        return self.batch_insert_function(to_insert)


    @staticmethod
    def _resolve_sources(source, n):
        """Normalizes `source` (a single int, or a list/tuple of `n` ints) into a list of `n` ints.

        Args:
            source: an integer provenance index applied to every entry, or a
                    list/tuple of `n` integer provenance indices, one per entry.
            n: the number of entries being inserted.

        Returns:
            A list of `n` integers.

        Raises:
            ValueError: if `source` is a list/tuple whose length differs from `n`.
        """
        if isinstance(source, (list, tuple)):
            if len(source) != n:
                raise ValueError(
                    "expected a list of {} sources (one per function), got {}".format(
                        n, len(source)
                    )
                )
            return list(source)
        return [source] * n


    def insert_quadratic(self, s, source=-1, thk_spec=None, ccz_size=-1, aut_size=-1, comment=None):
        """Inserts a quadratic APN function as a TYPE_QUADRATIC entry in a new CCZ class.

        Args:
            s: an S-boxable quadratic APN function.
            source: integer index into biblio_dict.from_which_paper(self.n, source).
            thk_spec: pre-computed thickness spectrum of `s` (e.g. from
                `get_WalshZeroesSpaces(s).thickness_spectrum()`, already needed by
                callers such as `add_quadratic_ccz_class` to transport thicknesses
                to the rest of the CCZ class), or None to compute it here directly.
            ccz_size: number of EA classes in this function's CCZ class, if already
                known (e.g. by `add_quadratic_ccz_class`), or -1 if unknown.
            aut_size: size of Aut(s), if already known, or -1 if unknown.
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            The row id of the inserted entry.
        """
        sb = get_sbox(s)
        if algebraic_degree(sb) != 2:
            raise ValueError("insert_quadratic: function is not quadratic")
        mug = apn_ea_mugshot(sb)
        invariant = self._make_invariant(mug)
        rep = bytearray(quadratic_compact_representation(sb.lut()))
        lin = linearity(sb)
        if thk_spec is None:
            thk_spec = thickness_spectrum(sb)
        thk = thk_spec.maximum()
        inserted_id = self._insert_entry(
            TYPE_QUADRATIC, self.number_of_ccz_classes, invariant, rep,
            2, lin, thk, source, ccz_size=ccz_size, aut_size=aut_size
        )
        self.number_of_ccz_classes += 1
        self.log_journal("insert_quadratic", comment)
        return inserted_id


    def insert_ccz_mapping(self, quad_ccz_id, mapping, s, source=-1, thk=None, mug=None, comment=None):
        """Inserts a CCZ-equivalent-to-quadratic function as a TYPE_CCZ_MAPPING entry.

        Args:
            quad_ccz_id: the ccz_id of the TYPE_QUADRATIC reference entry.
            mapping: the F2AffineMap L such that ccz_equivalent_function(quad_ref, L) == s.
            s: the S-boxable function to store (used only to compute the invariant).
            source: integer provenance index.
            thk: pre-computed thickness (see `populate_quadratic_ccz_class`, which derives
                 it cheaply from the reference's Walsh zero spaces instead of recomputing
                 them from scratch for `s`), or None to compute it here the direct way.
            mug: pre-computed EA mugshot (see `populate_quadratic_ccz_class`, which builds it
                 from spectra partly reused from the quadratic reference instead of recomputing
                 everything from `s`), or None to compute it here the direct way.
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            The row id of the inserted entry.
        """
        sb = get_sbox(s)
        if mug is None:
            mug = apn_ea_mugshot(sb)
        invariant = self._make_invariant(mug)
        rep = self._mapping_to_bytes(mapping)
        deg = algebraic_degree(sb)
        lin = linearity(sb)
        if thk is None:
            thk = thickness_spectrum(sb).maximum()
        inserted_id = self._insert_entry(
            TYPE_CCZ_MAPPING, int(quad_ccz_id), invariant, rep,
            deg, lin, thk, source
        )
        self.log_journal("insert_ccz_mapping", comment)
        return inserted_id


    def insert_non_quadratic(self, s, source=-1, comment=None):
        """Inserts a non-CCZ-quadratic APN function as a TYPE_NONQUADRATIC entry in a new CCZ class.

        Args:
            s: an S-boxable APN function that is not CCZ-equivalent to any quadratic.
            source: integer provenance index.
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            The row id of the inserted entry.
        """
        sb = get_sbox(s)
        abs_walsh = absolute_walsh_spectrum(sb)
        deg_spec   = degree_spectrum(sb)
        thk_spec   = thickness_spectrum(sb)
        mug = apn_ea_mugshot_from_spectra(
            abs_walsh,
            deg_spec,
            sigma_multiplicities(sb, k=4),
            thk_spec
        )
        invariant = self._make_invariant(mug)
        rep = bytearray(sb.to_bytes())
        deg = deg_spec.maximum()
        lin = abs_walsh.maximum()
        thk = thk_spec.maximum()
        inserted_id = self._insert_entry(
            TYPE_NONQUADRATIC, self.number_of_ccz_classes, invariant, rep,
            deg, lin, thk, source
        )
        self.number_of_ccz_classes += 1
        self.log_journal("insert_non_quadratic", comment)
        return inserted_id


    def batch_insert_quadratic_function(self, functions, source=-1, comment=None):
        """Inserts many quadratic APN functions in a single bulk operation, each
        as a TYPE_QUADRATIC entry in its own new CCZ class.

        Equivalent to calling `insert_quadratic` once per function, but performs
        a single bulk INSERT instead of one INSERT per function.

        The "thickness" column is not computed here and is stored as
        `_THICKNESS_NOT_COMPUTED` (-1): computing it requires enumerating all
        Walsh zero spaces of the function for efficiency reasons. Use the populate function
        or the update_database function to add thel afterwards.

        Args:
            functions: an iterable of S-boxable quadratic APN functions.
            source: either a single integer index into
                    biblio_dict.from_which_paper(self.n, source), applied to every
                    inserted entry, or a list/tuple of such indices, one per
                    element of `functions` (must have the same length).
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            A list of row ids, one per element of `functions`, in the same order.

        Raises:
            ValueError: if `source` is a list/tuple whose length differs from
                        the number of functions.
        """
        functions = list(functions)
        sources = self._resolve_sources(source, len(functions))
        ccz_id0 = self.number_of_ccz_classes
        rows = []
        for i, s in enumerate(functions):
            sb = get_sbox(s)
            if algebraic_degree(sb) != 2:
                raise ValueError(
                    "batch_insert_quadratic_function: function at index {} is not quadratic".format(i)
                )
            mug = apn_ea_mugshot(sb)
            invariant = self._make_invariant(mug)
            rep = quadratic_compact_representation(sb.lut())
            rows.append((
                TYPE_QUADRATIC, ccz_id0 + i, invariant, rep,
                2, linearity(sb), self._THICKNESS_NOT_COMPUTED, sources[i], -1, -1
            ))
        ids = self._insert_entries(rows)
        self.number_of_ccz_classes += len(rows)
        self.log_journal("batch_insert_quadratic_function", comment)
        return ids


    def batch_insert_non_quadratic_function(self, functions, source=-1, comment=None):
        """Inserts many non-CCZ-quadratic APN functions in a single bulk operation,
        each as a TYPE_NONQUADRATIC entry in its own new CCZ class.

        Equivalent to calling `insert_non_quadratic` once per function, but
        performs a single bulk INSERT instead of one INSERT per function.

        Args:
            functions: an iterable of S-boxable APN functions, none of which is
                       CCZ-equivalent to a quadratic.
            source: either a single integer provenance index, applied to every
                    inserted entry, or a list/tuple of such indices, one per
                    element of `functions` (must have the same length).
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            A list of row ids, one per element of `functions`, in the same order.

        Raises:
            ValueError: if `source` is a list/tuple whose length differs from
                        the number of functions.
        """
        functions = list(functions)
        sources = self._resolve_sources(source, len(functions))
        ccz_id0 = self.number_of_ccz_classes
        rows = []
        for i, s in enumerate(functions):
            sb = get_sbox(s)
            abs_walsh = absolute_walsh_spectrum(sb)
            deg_spec  = degree_spectrum(sb)
            thk_spec  = thickness_spectrum(sb)
            mug = apn_ea_mugshot_from_spectra(
                abs_walsh, deg_spec, sigma_multiplicities(sb, k=4), thk_spec
            )
            invariant = self._make_invariant(mug)
            rows.append((
                TYPE_NONQUADRATIC, ccz_id0 + i, invariant, sb.to_bytes(),
                deg_spec.maximum(), abs_walsh.maximum(), thk_spec.maximum(), sources[i], -1, -1
            ))
        ids = self._insert_entries(rows)
        self.number_of_ccz_classes += len(rows)
        self.log_journal("batch_insert_non_quadratic_function", comment)
        return ids


    def update_database(self, indices, field, values, comment=None):
        """Updates a single field for a list of entries.

        Args:
            indices: list of row ids to update.
            field:   the column name to set (must be a key of self.row_structure).
            values:  list of new values, parallel to `indices`.
            comment: free-text note describing why this update is being made, logged
                      to the journal (see `log_journal`), or None.

        Raises:
            ValueError: if `field` is not a column of this database, or if the
                        lengths of `indices` and `values` differ.
        """
        if field not in self.row_structure:
            raise ValueError(
                "update_database: unknown field '{}' (valid fields: {})".format(
                    field, list(self.row_structure.keys())
                )
            )
        if len(indices) != len(values):
            raise ValueError(
                "update_database: indices and values must have the same length "
                "({} vs {})".format(len(indices), len(values))
            )
        query = "UPDATE {} SET {} = ? WHERE id = ?".format(self.functions_table, field)
        for idx, val in zip(indices, values):
            self.cursor.execute(query, (val, int(idx)))
        self.log_journal("update_database", comment)


    def populate_quadratic_ccz_class(self, entry_id, comment=None):
        """Populates an existing TYPE_QUADRATIC entry with its CCZ-equivalent EA classes.

        Finds the TYPE_QUADRATIC entry with the given `entry_id`, computes the
        Walsh-zero-space mapping for its CCZ class, and inserts each EA-class
        representative as a TYPE_CCZ_MAPPING entry — but only if the class has not
        already been populated it is checked by verifying if the number of functions
        with the same ccz_id is greater than 1. It is a partial test, so always verify
        the state of the database before calling this function.

        The source of every inserted mapping entry is inherited from the TYPE_QUADRATIC entry.

        Args:
            entry_id: the row id of a TYPE_QUADRATIC entry.
            comment: free-text note logged to the journal (see `log_journal`) for every
                     modification made by this call (the thickness/ccz_size/aut_size field
                     updates and the batch insert of EA representatives), or None.

        Returns:
            A list of the inserted row ids, or an empty list if the class was already
            fully populated.

        Raises:
            ValueError: if no entry with the given id exists, or if it is not TYPE_QUADRATIC.
        """
        entries = self.query_functions({"id": entry_id})
        if not entries:
            raise ValueError(
                "populate_quadratic_ccz_class: no entry with id={}".format(entry_id)
            )
        entry = entries[0]
        if entry["type"] != TYPE_QUADRATIC:
            raise ValueError(
                "populate_quadratic_ccz_class: entry {} has type {} (expected TYPE_QUADRATIC={})".format(
                    entry_id, entry["type"], TYPE_QUADRATIC
                )
            )

        sb = entry["sbox"]
        ccz_id = entry["ccz_id"]
        source = entry["source"]

        n_existing = len(self.query_functions({"ccz_id": ccz_id}))

        # Populating only if there is one function
        if n_existing > 1:
            return []

        # Computing Walsh Zeroes and adding the thickness to the quadratic
        ws_full = get_WalshZeroesSpaces(sb)
        self.update_database([entry_id], "thickness", [ws_full.thickness_spectrum().maximum()], comment= None)

        # Absolute Walsh Spectrum is CCZ-invariant, we compmute it once
        abs_walsh = absolute_walsh_spectrum(sb)
        # linearity is also CCZ-invariant
        lin = abs_walsh.maximum()

        # Aut(sb) and the orbit reduction it induces are computed explicitly, instead of
        # calling get_WalshZeroesSpaces_quadratic_apn, which would redo the Walsh zero space
        # search from scratch
        aut = automorphisms_from_ortho_derivative(sb)
        ws_reduced = ws_full.copy()
        ws_reduced.init_mappings_using_automorphisms(aut)
        mappings = ws_reduced.get_mappings()
        aut_size = len(aut)
        ccz_size = len(mappings)
        self.update_database([entry_id], "ccz_size", [ccz_size], comment=None)
        self.update_database([entry_id], "aut_size", [aut_size], comment=None)

        # Rows are collected and inserted with a single executemany (via _insert_entries)
        rows = []
        valve = False
        for L in mappings:
            g = ccz_equivalent_function(sb, L)
            if valve == False:
                if algebraic_degree(g) == 2:
                    # quadratic results are covered by the TYPE_QUADRATIC entry so we continue
                    # and we change valve to True so we do not compute the degree afterwards
                    valve = True
                    continue
            # We transport the Walsh Zeroes as computing them from scratch is the costliest
            thk_spec = ws_full.image_by(L.inverse().transpose()).thickness_spectrum()
            # We compute the rest of the mugshot
            deg_spec = degree_spectrum(g)
            mug = apn_ea_mugshot_from_spectra(
                abs_walsh, deg_spec, sigma_multiplicities(g, k=4), thk_spec
            )
            rows.append((
                TYPE_CCZ_MAPPING, ccz_id, self._make_invariant(mug), self._mapping_to_bytes(L),
                deg_spec.maximum(), lin, thk_spec.maximum(), source, ccz_size, aut_size
            ))
        inserted_ids = self._insert_entries(rows)
        self.log_journal("populate_quadratic_ccz_class", comment)
        return inserted_ids


    def populate_non_quadratic_ccz_class(self, entry_id, comment=None):
        """Populates an existing TYPE_NONQUADRATIC entry with its CCZ-class EA representatives.

        The reference's graph automorphism group Aut(f) is computed via `ccz_equivalences(f, f)` ,
        a generic LAT-based self-equivalence search.

        Population is skipped when the number of entries sharing the reference's ccz_id is
        already greater than 1, meaning the class has been (at least partially) populated
        before. It is a partial test, so always verify
        the state of the database before calling this function.

        The source of every inserted entry is inherited from the TYPE_NONQUADRATIC reference.

        Args:
            entry_id: the row id of a TYPE_NONQUADRATIC entry.
            comment: free-text note logged to the journal (see `log_journal`) for every
                     modification made by this call (the ccz_size/aut_size field updates and
                     every inserted EA representative), or None.

        Returns:
            A list of the inserted row ids, or an empty list if the class was already
            (at least partially) populated.

        Raises:
            ValueError: if no entry with the given id exists, or if it is not TYPE_NONQUADRATIC.
        """
        entries = self.query_functions({"id": entry_id})
        if not entries:
            raise ValueError(
                "populate_non_quadratic_ccz_class: no entry with id={}".format(entry_id)
            )
        entry = entries[0]
        if entry["type"] != TYPE_NONQUADRATIC:
            raise ValueError(
                "populate_non_quadratic_ccz_class: entry {} has type {} "
                "(expected TYPE_NONQUADRATIC={})".format(
                    entry_id, entry["type"], TYPE_NONQUADRATIC
                )
            )

        sb = entry["sbox"]
        ccz_id = entry["ccz_id"]
        source = entry["source"]

        n_existing = len(self.query_functions({"ccz_id": ccz_id}))
        if n_existing > 1:
            return []

        # Graph automorphisms via a generic self-equivalence search
        aut_pairs = ccz_equivalences(sb, sb)
        aut_size = len(aut_pairs)
        aut_linear_parts = [A for A, cst in aut_pairs]

        ws_reduced = get_WalshZeroesSpaces(sb).copy()
        ws_reduced.init_mappings_using_automorphisms(aut_linear_parts)
        mappings = ws_reduced.get_mappings()
        ccz_size = len(mappings)

        self.update_database([entry_id], "ccz_size", [ccz_size], comment=None)
        self.update_database([entry_id], "aut_size", [aut_size], comment=None)

        inserted_ids = []
        for L in mappings:
            g_sb = ccz_equivalent_function(sb, L)
            if len(g_sb) == 0:
                continue
            _, mug = self._compute_invariant(g_sb)
            if not self.is_new(g_sb, mug=mug):
                continue
            invariant = bytearray(self._make_invariant(mug))
            rep = bytearray(g_sb.to_bytes())
            deg = algebraic_degree(g_sb)
            lin = linearity(g_sb)
            thk = thickness_spectrum(g_sb).maximum()
            inserted_ids.append(
                self._insert_entry(
                    TYPE_NONQUADRATIC, ccz_id, invariant, rep, deg, lin, thk, source,
                    ccz_size=ccz_size, aut_size=aut_size
                )
            )
        self.log_journal("populate_non_quadratic_ccz_class", comment)
        return inserted_ids


    def add_quadratic_ccz_class(self, s, source=-1, comment=None):
        """Inserts a complete CCZ-equivalence class of a quadratic APN function.

        Args:
            s: an S-boxable quadratic APN function.
            source: integer provenance index (applied to all inserted entries).
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            A list of row ids: [quad_id, mapping_id_0, mapping_id_1, ...].
        """
        sb = get_sbox(s)
        if algebraic_degree(sb) != 2:
            raise ValueError("add_quadratic_ccz_class: function is not quadratic")

        # !CAREFUL! We use the normalized form composed ONLY of the quadratic ANF terms
        sb_stored = get_sbox(
            quadratic_sbox_from_compact_representation(
                quadratic_compact_representation(sb.lut()),
                sb.get_input_length(),
                sb.get_output_length()
            )
        )

        # Save the ccz_id that insert_quadratic will use before it increments the counter.
        quad_ccz_id = self.number_of_ccz_classes

        # Computing Walsh Zeroes once and reusing them (via transport) for every EA class
        ws_full = get_WalshZeroesSpaces(sb_stored)

        # Aut(sb_stored) and the orbit reduction it induces are computed explicitly,instead
        # of calling get_WalshZeroesSpaces_quadratic_apn, which would redo the Walsh zero
        # space search from scratch
        aut = automorphisms_from_ortho_derivative(sb_stored)
        ws_reduced = ws_full.copy()
        ws_reduced.init_mappings_using_automorphisms(aut)
        mappings = ws_reduced.get_mappings()
        aut_size = len(aut)
        ccz_size = len(mappings)

        quad_id = self.insert_quadratic(
            sb_stored, source, thk_spec=ws_full.thickness_spectrum(),
            ccz_size=ccz_size, aut_size=aut_size, comment=comment
        )
        inserted_ids = [quad_id]

        # Absolute Walsh Spectrum is CCZ-invariant, we compute it once
        abs_walsh = absolute_walsh_spectrum(sb_stored)
        # linearity is CCZ-invariant
        lin = abs_walsh.maximum()

        # Rows are collected and inserted with a single executemany (via _insert_entries)
        rows = []
        valve = False
        for L in mappings:
            g = ccz_equivalent_function(sb_stored, L)
            if valve == False:
                if algebraic_degree(g) == 2:
                    # quadratic results are covered by the TYPE_QUADRATIC entry so we continue
                    # and we change valve to True so we do not compute the degree afterwards
                    valve = True
                    continue
            # We transport the Walsh Zeroes as computing them from scratch is the costliest
            thk_spec = ws_full.image_by(L.inverse().transpose()).thickness_spectrum()
            # deg_spec.maximum() gives the degree without a separate algebraic_degree(g)
            deg_spec = degree_spectrum(g)
            mug = apn_ea_mugshot_from_spectra(
                abs_walsh, deg_spec, sigma_multiplicities(g, k=4), thk_spec
            )
            rows.append((
                TYPE_CCZ_MAPPING, quad_ccz_id, self._make_invariant(mug), self._mapping_to_bytes(L),
                deg_spec.maximum(), lin, thk_spec.maximum(), source, ccz_size, aut_size
            ))
        inserted_ids += self._insert_entries(rows)
        self.log_journal("add_quadratic_ccz_class", comment)
        return inserted_ids


    def add_non_quadratic_ccz_class(self, s, source=-1, comment=None):
        """Inserts a complete CCZ-equivalence class rooted at a non-quadratic APN function.

        Inserts `s` as a TYPE_NONQUADRATIC reference entry, then enumerates all admissible
        mappings without filtering (automorphisms are not available for non-quadratic functions).
        Each CCZ-equivalent function that is EA-new with respect to the current database state
        is inserted as a TYPE_NONQUADRATIC entry under the same ccz_id.  EA-newness is checked
        via `is_new` after every insertion so that already-inserted representatives are not
        duplicated by later mappings in the same run.

        Args:
            s: an S-boxable APN function of algebraic degree strictly greater than 2.
            source: integer index into `from_which_paper(self.n, source)` for provenance.
            comment: free-text note logged to the journal (see `log_journal`), or None.

        Returns:
            A list of row ids: [reference_id, new_id_0, new_id_1, ...].

        Raises:
            ValueError: if `s` is quadratic (use `add_quadratic_ccz_class` instead).
        """
        sb = get_sbox(s)
        if algebraic_degree(sb) == 2:
            raise ValueError(
                "add_non_quadratic_ccz_class: function is quadratic; "
                "use add_quadratic_ccz_class instead"
            )

        new_ccz_id = self.number_of_ccz_classes
        ref_id = self.insert_non_quadratic(sb, source, comment=comment)

        ws = get_WalshZeroesSpaces(sb)
        inserted_ids = [ref_id]
        for L in ws.get_mappings():
            g_sb = ccz_equivalent_function(sb, L)
            if len(g_sb) == 0:
                continue
            _, mug = self._compute_invariant(g_sb)
            if not self.is_new(g_sb, mug=mug):
                continue
            invariant = bytearray(self._make_invariant(mug))
            rep = bytearray(g_sb.to_bytes())
            deg = algebraic_degree(g_sb)
            lin = linearity(g_sb)
            thk = thickness_spectrum(g_sb).maximum()
            inserted_ids.append(
                self._insert_entry(TYPE_NONQUADRATIC, new_ccz_id, invariant, rep, deg, lin, thk, source)
            )
        self.log_journal("add_non_quadratic_ccz_class", comment)
        return inserted_ids
