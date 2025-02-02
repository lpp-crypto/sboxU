#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont
import platform

from .sboxU_cython import *
from .utils import *
from .display import *
from .diff_lin import *
from .ccz import *
from .cycles import *


# !SECTION! Statistical anomalies 
# -------------------------------

class StatisticalAnomaly:
    """A class that is not supposed to be used directly, but from
    which several anomalies inherit.

    """
    def __init__(self, s, name, spectrum):
        self.name = name
        self.n = int(lg2(len(s)))
        self.spectrum = spectrum


    def summary(self):
        return [
            "spectrum:   {}".format(pretty_spectrum(self.spectrum)),
            "positive A: {:10.4f}".format(float(self.positive_anomaly)),
            "negative A: {:10.4f}".format(float(self.negative_anomaly)),
        ]

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"
        return result


class DifferentialAnomaly(StatisticalAnomaly):
    """The differential anomaly quantifies the probability that a
    random function has a pair (u, n(u)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, u is the differential uniformity, and n(u) is
    its number of occurrences in the DDT.

    """
    def __init__(self, s):
        StatisticalAnomaly.__init__(self, s, "differential", differential_spectrum(s))
        if max(self.spectrum.keys()) == 2:
            # APN-ness basically doesn't happen
            self.positive_anomaly = -infinity
            self.negative_anomaly = +infinity
        else:
            self.positive_anomaly = table_anomaly(s, "DDT", spec=self.spectrum)
            self.negative_anomaly = table_negative_anomaly(s, "DDT", spec=self.spectrum)
        
    
class LinearAnomaly(StatisticalAnomaly):
    """The linear anomaly quantifies the probability that a
    random function has a pair (l, n(l)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, l is the linearity, and n(l) is
    the sum of the number of occurrences of l and -l in the LAT.

    """
    def __init__(self, s):
        spec = walsh_spectrum(s)
        abs_spec = defaultdict(int)
        for k in spec.keys():
            abs_spec[abs(k)] += spec[k]
        StatisticalAnomaly.__init__(self, s, "linear", dict(abs_spec))
        self.positive_anomaly = table_anomaly(s, "LAT", spec=self.spectrum)
        self.negative_anomaly = table_negative_anomaly(s, "LAT", spec=self.spectrum)

        
class BoomerangAnomaly(StatisticalAnomaly):
    """The boomerang anomaly quantifies the probability that a
    random function has a pair (b, n(b)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, b is the boomerang uniformity, and n(b) is
    its number of occurrences in the BCT.

    """
    def __init__(self, s):
        StatisticalAnomaly.__init__(self, s, "boomerang", boomerang_spectrum(s))
        self.positive_anomaly = table_anomaly(s, "BCT", spec=self.spectrum)
        self.negative_anomaly = table_negative_anomaly(s, "BCT", spec=self.spectrum)


def card_vectors_rank(r, l, n):
    """Returns the number of families of vectors of `n` bits
    containing `l` vectors and of rank `r`.

    """
    if r > l:
        # rank can't be higher than the number of vectors
        return 0 
    elif r == 0:
        # exactly 1 family of rank 0: the all-zero one
        return 1 
    else:
        # a family of rank r and length l can be built in two distinct
        # ways:
        # - from a family of rank r-1 and length l-1 where we add a
        #   vector not in their span, or
        # - from a family of rank r and length l-1 where we add a
        #   vector in their span.
        return (2**n - 2**(r-1)) * card_vectors_rank(r-1, l-1, n) + 2**r * card_vectors_rank(r, l-1, n) 


class DegreeAnomaly(StatisticalAnomaly):
    """The degree anomaly quantifies the probability that a
    random function has a pair (d, n(d)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, d is the maximum degree of the components of
    the function, and n(b) is the number of components with this
    degree.

    """
    def __init__(self, s):
        StatisticalAnomaly.__init__(self, s, "degree", degree_spectrum(s))
        self.is_perm = is_permutation(s)
        self.hdim = hdim(s)
        self.hdim_rank = self.hdim.rank()
        ranks_card = [card_vectors_rank(r, self.n, self.n) for r in range(0, self.n+1)]
        card_worse = sum(ranks_card[r] for r in range(0, self.hdim_rank))
        if card_worse > 0:
            self.negative_anomaly = self.n**2 - lg2(card_worse)
        else:
            self.negative_anomaly = self.n * 2**self.n
        card_better = sum(ranks_card[r] for r in range(self.hdim_rank, self.n+1))
        if card_better > 0:
            self.positive_anomaly = self.n**2 - lg2(card_better)
        else:
            self.positive_anomaly = self.n * 2**self.n
        if not is_permutation(s):
            self.negative_anomaly += self.spectrum[self.n]

    def summary(self):
        result = StatisticalAnomaly.summary(self)
        result.append("HDIM: [{}]".format(str(self.hdim[0])[1:-1]))
        for line in self.hdim[1:]:
            result.append("      [{}]".format(str(line)[1:-1]))
        result.append("rank(HDIM) = {}".format(self.hdim_rank))
        return result

    
# !SECTION! Structural anomalies
# ------------------------------

class CPSAnomaly:
    def __init__(self, s):
        self.name = "CPS"
        preimages = defaultdict(int)
        for y in s:
            preimages[y] += 1
        self.spectrum = defaultdict(int)
        for y in preimages.keys():
            self.spectrum[preimages[y]] += 1
        self.positive_anomaly = 0
        self.negative_anomaly = 0

    def summary(self):
        return [
            "CPS : {}".format(pretty_spectrum(self.spectrum)),
        ]

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"


class LinearStructAnomaly:
    def __init__(self, s):
        self.name = "linear structures"
        self.struct = linear_structures_vectorial(s)
        self.spectrum = defaultdict(int)
        for c in self.struct.keys():
            self.spectrum[len(self.struct[c][0]) + len(self.struct[c][1])] += 1
        self.positive_anomaly = 0
        self.negative_anomaly = 0

        
    def summary(self):
        result = [
            "Linear Structures : {}".format(pretty_spectrum(self.spectrum)),
        ]
        for c in sorted(self.struct.keys()):
            result.append(" - {:2x} : 0:{}".format(c, pretty_vector(self.struct[c][0])))
            result.append("        1:{}".format(pretty_vector(self.struct[c][1])))
        return result
    

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"

            
class CycleAnomaly:
    def __init__(self, s):
        self.name = "cycle decomposition"
        self.perm = is_permutation(s)
        if self.perm:
            self.cycles = cycle_decomposition(s)
            self.spectrum = defaultdict(int)
            for c in self.cycles:
                self.spectrum[len(c)] += 1
        else:
            self.spectrum = {}
        # !TODO! evaluate the anomaly of all cycle types 
        self.positive_anomaly = 0
        self.negative_anomaly = 0

    def summary(self):
        if not self.perm:
            return ["Not a permutation"]
        result = ["Cycle type : {}".format(pretty_spectrum(self.spectrum))]
        cycles_by_length = defaultdict(list)
        for c in self.cycles:
            cycles_by_length[len(c)].append(c)
        # cycles of a specific length (1 or 2)
        if len(cycles_by_length[1]) > 1:
            fixed_points_str = "fixed points   (#={:d}) : [".format(len(cycles_by_length[1]))
            for c in cycles_by_length[1]:
                fixed_points_str += "{}, ".format(c[0])
            result.append(fixed_points_str[:-2] + "]")
        else:
            result.append("fixed points   (#=0) : []")
        if len(cycles_by_length[2]) > 0:
            transpositions_str = "transpositions (#={:d}) : [".format(len(cycles_by_length[2]))
            for c in cycles_by_length[2]:
                transpositions_str += "{}, ".format(c)
            result.append(transpositions_str[:-2] + "]")
        else:
            result.append("transpositions (#=0) : []")
        # generic cycles
        for l in sorted(cycles_by_length.keys()):
            if l > 2:
                for c in cycles_by_length[l]:
                    result.append("- [len={:3d}]  {}".format(len(c), c))
        return result

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"
    
    

class CCZAnomaly:
    def __init__(self, s, threshold=20, spaces=None):
        self.name = "CCZ"
        N = int(log(len(s), 2))
        self.lat_zeroes = lat_zeroes(s)
        self.spaces_all_found = True
        if spaces == None:
            self.spaces = []
            # doing the thickness spectrum by hand so that this does not
            # take up too much time a priori
            for space in vector_spaces_bases_iterator(self.lat_zeroes, N, 2*N):
                self.spaces.append(space)
                if len(self.spaces) == threshold:
                    self.spaces_all_found = False
                    break
        else:
            self.spaces = spaces
        self.spectrum = thickness_spectrum(s, spaces=self.spaces)
        self.negative_anomaly = 0
        self.positive_anomaly = 0
        for t in self.spectrum.keys():
            if t > 0:
                log_tu_density = N*2**N - 2**(N-t) * sum(lg2(i) for i in range(1, 2**t+1)) - t * (N-t) * 2**(N-t)
                log_affine_card = N + sum(lg2(2**N - 2**i) for i in range(0, N))
                multiplier = self.spectrum[t]
                if t == N and is_permutation(s):
                    # in a permutation, having one space of thickness n is normal
                    multiplier -= 1
                self.positive_anomaly += (log_tu_density - 2*log_affine_card - N**2) * multiplier


    def summary(self):
        result = []
        if not self.spaces_all_found:
            result += [
                "/!\\ WARNING: not all spaces found!",
                "    Rerun with a higher threshold to have them all."
            ]
        result += [
            "|Walsh zero spaces| = {}".format(len(self.spaces)),
            "thickness spectrum = {}".format(pretty_spectrum(self.spectrum)),
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
        ]
        return result
    

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"



# !SECTION! Automated analysis 
# ----------------------------

def affine_equivalence_monomial(s, d):
    """Assuming `s` is the LUT of a function and `d` is its DDT,
    checks if `s` could be affine equivalent to a monomial.

    """
    first_row_count = defaultdict(int)
    for b in d[1]:
        first_row_count[b] += 1
    for a in range(2, len(d)):
        ath_row_count = defaultdict(int)
        for b in d[a]:
            ath_row_count[b] += 1
            if ath_row_count[b] > first_row_count[b]:
                return False
    return True
        

def sort_table_rows(t, threshold=2):
    row_spectra = defaultdict(list)
    for a in range(1, len(t)):
        row = t[a]
        spec = defaultdict(int)
        for b in range(1, len(row)):
            spec[abs(row[b])] += 1
        row_spectra[pretty_spectrum(spec)].append(a)
    result = {}
    for spec_str in sorted(row_spectra.keys()):
        if len(row_spectra[spec_str]) > threshold:
            result[spec_str] = row_spectra[spec_str]
    return result
    

# !SUBSECTION! Main class
# -----------------------

ALWAYS_ANALYSED_SIZE = 14

class Analysis:
    def __init__(self,
                 s,
                 name="sbox",
                 # which analyses to perform
                 differential=None,
                 linear=None,
                 boomerang=None,
                 degree=None,
                 cps=None,
                 cycles=None,
                 linear_struct=None,
                 ccz=None,
                 tu=None,
                 # how to perform the analysis
                 store_tables=None,
                 deep=False,
                 textures=False,
            ):
        if isinstance(s, sage.crypto.sbox.SBox):
            self.N = s.n
            s = list(s)
        elif isinstance(s, list):
            self.N = int(lg2(len(s)))
        else:
            raise Exception("unsupported input type: {}".format(type(s)))
        self.name = name
        self.deep = deep
        if isinstance(s, sage.crypto.sbox.SBox):
            self.lut = list(s)
        elif isinstance(s, list):
            self.lut = s
        else:
            raise Exception("s must be a list (or an SBox instance)")
        # setting up parameters of the analysis
        if cps == None:
            cps = True
        if cycles == None:
            cycles = True
        if differential == None:
            differential = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if linear == None:
            linear = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if boomerang == None:
            boomerang = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if degree == None:
            degree = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if linear_struct == None:
            linear_struct = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if ccz == None:
            ccz = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        if store_tables == None:
            store_tables = (self.N <= ALWAYS_ANALYSED_SIZE) or deep
        self.textures = textures
        self.ccz = ccz
        # handling basic properties
        self.anomalies = {}
        self.tables = {}
        self.spectra = {}
        if differential:
            self.anomalies["DDT"] = DifferentialAnomaly(s)
            self.spectra["DDT"] = self.anomalies["DDT"].spectrum
            if store_tables:
                self.tables["DDT"] = ddt(s)
        if linear:
            self.anomalies["LAT"] = LinearAnomaly(s)
            self.spectra["LAT"] = self.anomalies["LAT"].spectrum
            if store_tables:
                self.tables["LAT"] = lat(s)
        if boomerang:
            self.anomalies["BCT"] = BoomerangAnomaly(s)
            self.spectra["BCT"] = self.anomalies["BCT"].spectrum
            if store_tables:
                self.tables["BCT"] = bct(s)
        if degree:
            self.anomalies["DEG"] = DegreeAnomaly(s)
            self.spectra["DEG"] = self.anomalies["DEG"].spectrum
        if cps:
            self.anomalies["CPS"] = CPSAnomaly(s)
            self.spectra["CPS"] = self.anomalies["CPS"].spectrum
        if linear_struct:
            self.anomalies["LST"] = LinearStructAnomaly(s)
            self.spectra["LST"] = self.anomalies["LST"].spectrum
        if ccz:
            if deep:
                self.anomalies["CCZ"] = CCZAnomaly(s, threshold=Infinity)
            else:
                self.anomalies["CCZ"] = CCZAnomaly(s)
            self.spectra["CCZ"] = self.anomalies["CCZ"].spectrum
        if cycles:
            self.anomalies["CYC"] = CycleAnomaly(s)
        if tu:
            # !CONTINUE!  Finish automated TU analysis
            print("Unsupported yet")
        self.advanced = {}
        if deep:
            self.advanced_analysis()


    def advanced_analysis(self):
        # more sophisticated analysis
        # !TODO! advanced analysis
        # Looking again at the tables: identical lines/rows and textures
        if self.textures:
            self.advanced["textures"] = {}
        self.advanced["common spectra"] = {table_name : defaultdict(list)
                                           for table_name in self.tables.keys()}
        for table_name in sorted(self.tables.keys()):
            t = self.tables[table_name]
            threshold = 8 if table_name == "DDT" else 2
            row_spectra = sort_table_rows(t, threshold=threshold)
            col_spectra = sort_table_rows([[t[a][b] for a in range(0, 2**self.N)]
                                           for b in range(0, 2**self.N)],
                                          threshold=threshold)
            if len(row_spectra.keys()) > 0:
                self.advanced["common spectra"][table_name]["rows"] = row_spectra
            if len(col_spectra.keys()) > 0:
                self.advanced["common spectra"][table_name]["cols"] = col_spectra
            if self.textures:
                self.advanced["textures"][table_name] = {
                    "ADD" : add_texture(self.tables[table_name]),
                    "XOR" : xor_texture(self.tables[table_name])
                }
        # CCZ properties: TU-decompositions, etc.
        if not self.ccz:
            # if ccz was not set then we need to look for all the
            # vector spaces of zeroes
            self.anomalies["CCZ"] = CCZAnomaly(s, threshold=Infinity, spaces=get_lat_zeroes(self.lut))
            self.spectra["CCZ"] = self.anomalies["CCZ"].spectrum
        if self.anomalies["CCZ"].spectrum[self.N] > 0:
            self.advanced["EA permutations"] = ea_equivalent_permutation_mappings(
                self.lut, 
                self.anomalies["CCZ"].spaces
            )
          
            
    def expected_distribution(self, table_name):
        return get_proba_func(self.lut, table_name)

    
    def show(self,
             indent="",    
             noteworthy_threshold=10,
             noteworthy_mark="  # Noteworty"):
        print(indent + "LUT:")
        print(indent + "      {}".format(
            pretty_vector(list(range(0, min(2**self.N, 16))), template="{:02x} ")[1:-1]
        ))
        print(indent + "-"*54)
        if self.N >= 4:
            for i in range(0, 2**(self.N-4)):
                print(indent + "{:2x}0 | {}".format(
                    i, pretty_vector(self.lut[i*16:i*16+16],
                                     template="{:02x} "
                                     )[1:-1]
                ))
        else:
            print(indent + " | {}".format(
                pretty_vector(self.lut[i*16:i*16+16],
                              template="{:02x} "
                              )[1:-1]
            ))
        good_properties = []
        bad_properties  = []
        for anomaly_name in sorted(self.anomalies.keys()):
            a = self.anomalies[anomaly_name]        
            to_print = a.name
            if a.positive_anomaly >= noteworthy_threshold:
                to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                good_properties.append(a.name)
            elif a.negative_anomaly >= noteworthy_threshold:
                to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                bad_properties.append(a.name)
            print("\n{}# {}\n".format(indent, to_print))
            for line in a.summary():
                print(indent + line)
        print("\n{}# {}\n".format(indent, "Preliminary Results:"))
        if len(good_properties) > 0:
            print(indent + "- Noteworthy qualities: " + str(good_properties)[1:-1])
        if len(bad_properties) > 0:
            print(indent + "- Noteworthy flaws    : " + str(bad_properties)[1:-1])
        if (len(good_properties) == 0) and (len(bad_properties) == 0):
            print(indent + "- The function is dull")
        if len(self.advanced.keys()) > 0:
            self.show_advanced(indent=indent)

            
    def show_advanced(self, indent=""):
        print("\n" + indent + "# Advanced properties")
        print("\n" + indent + "## Common spectra")
        # advanced table stuff
        for table_name in self.advanced["common spectra"]:
            if len(self.advanced["common spectra"][table_name]["rows"]) > 0:
                print("{}- {} rows : {}".format(
                    indent,
                    table_name,
                    self.advanced["common spectra"][table_name]["rows"]
                ))
            if len(self.advanced["common spectra"][table_name]["cols"]) > 0:
                print("{}- {} cols : {}".format(
                    indent,
                    table_name,
                    self.advanced["common spectra"][table_name]["cols"]
                ))
        # advanced CCZ
        if "EA permutations" in self.advanced.keys():
            print("\n" + indent + "## EA-equivalent permutations")
            for index, L in enumerate(self.advanced["EA permutations"]):
                print("\n- mapping number {:d}".format(index))
                for line in str(L).split("\n"):
                    print(indent + line)


    def save_pollock(self, 
                     name=None,
                     vmins=[None, None], 
                     vmaxs=[None, None], 
                     color_schemes=["CMRmap_r", "seismic"],
                     show=False,
                     cleanup=True):
        # pre-processing arguments
        if len(vmaxs) != len(color_schemes):
            raise Exception("Mismatched length in between `vmaxes` and `color_schemes`")
        elif len(vmaxs) != len(vmins):
            raise Exception("Mismatched length in between `vmaxes` and `vmins`")
        if name == None:
            name = self.name
        for x in [vmaxs, vmins, color_schemes]:
            if not isinstance(x, list):
                x = [x]
        # looping over all known properties
        for table_name in self.tables.keys():
            images_path = []
            spec = self.spectra[table_name]
            # Pollock representation(s)
            for i in range(0, len(vmaxs)):
                file_name = "{}-{}-{}".format(name, table_name, i)
                vmin = vmins[i] if vmins[i] != None else 0
                if vmaxs[i] != None:
                    vmax = vmaxs[i]  
                elif table_name == "DDT":
                    vmax = max(spec.keys())
                else:
                    avg_coeff = sum(spec[c]*abs(c) for c in spec.keys()) / sum(spec[c] for c in spec.keys() if c != 0)
                    vmax = min(floor(2*avg_coeff), max(spec.keys()))
                save_pollock(
                    self.tables[table_name],
                    name=file_name,
                    vmin=vmin,
                    vmax=vmax,
                    colorbar=True,
                    color_scheme=color_schemes[i],
                    title="Jackson Pollock Representation ({})".format(table_name),
                )
                images_path.append(file_name+".png")
            # textures
            if self.textures:
                file_name = "{}-xorTexture-{}".format(name, table_name)
                m = self.advanced["textures"][table_name]["XOR"]
                max_texture = max([m[i][j]
                                   for i, j in itertools.product(range(0, len(m)), range(0, len(m)))
                                   if (i,j) != (0, 0)])
                min_texture = min([m[i][j]
                                   for i, j in itertools.product(range(0, len(m)), range(0, len(m)))
                                   if (i,j) != (0, 0)])
                save_pollock(
                    m,
                    name=file_name,
                    vmin=min_texture,
                    vmax=max_texture,
                    colorbar=True,
                    color_scheme="CMRmap_r",
                    title="XOR texture".format(table_name),
                )
                images_path.append(file_name+".png")
                file_name = "{}-addTexture-{}".format(name, table_name)
                m = self.advanced["textures"][table_name]["ADD"]
                max_texture = max([m[i][j]
                                   for i, j in itertools.product(range(0, len(m)), range(0, len(m)))
                                   if (i,j) != (0, 0)])
                min_texture = min([m[i][j]
                                   for i, j in itertools.product(range(0, len(m)), range(0, len(m)))
                                   if (i,j) != (0, 0)])
                save_pollock(
                    m,
                    name=file_name,
                    vmin=min_texture,
                    vmax=max_texture,
                    colorbar=True,
                    color_scheme="CMRmap_r",
                    title="ADD texture".format(table_name),
                )
                images_path.append(file_name+".png")
            # comparison with expected distribution
            expected_distrib = self.expected_distribution(table_name)
            file_name = "distrib-{}-{}".format(name, table_name)
            plot_statistical(
                spec,
                n=self.N,
                expected_distrib=expected_distrib,
                file_name=file_name
            )
            images_path.append(file_name+".png")
            # patterns in the variance
            file_name = "{}-var-{}".format(name, table_name)
            plot_table_variances(self.tables[table_name], file_name=file_name)
            images_path.append(file_name+".png")
            # patterns by row and columns
            file_name = "distrib-by-row-{}-{}".format(name, table_name)
            plot_statistical_by_rows(
                self.tables[table_name],
                n=self.N,
                expected_distrib=self.expected_distribution(table_name),
                file_name=file_name
            )
            images_path.append(file_name+".png")
            file_name = "distrib-by-col-{}-{}".format(name, table_name)
            plot_statistical_by_rows(
                [[self.tables[table_name][a][b] for a in range(0, 2**self.N)]
                 for b in range(0, 2**self.N)],
                n=self.N,
                expected_distrib=self.expected_distribution(table_name) if table_name != "LAT" else lat_coeff_probability_permutation, # Léo: open problem: why is this necessary?
                file_name=file_name
            )
            images_path.append(file_name+".png")
            # assembling summary picture
            images = [Image.open(img_path) for img_path in images_path]
            width = max(sum(images[i].width for i in range(0, len(images), 2)),
                        sum(images[i].width for i in range(1, len(images), 2)))
            height_1 = max(images[i].height for i in range(0, len(images), 2))
            height_2 = max(images[i].height for i in range(1, len(images), 2))
            image_summary = Image.new("RGB", (width, height_1 + height_2))
            cursor_1 = 0
            cursor_2 = 0
            for i in range(0, len(images)):
                if (i % 2) == 0:
                    image_summary.paste(images[i], (cursor_1, 0))
                    cursor_1 += images[i].width
                else:
                    image_summary.paste(images[i], (cursor_2, height_1))
                    cursor_2 += images[i].width
            # writing stuff
            title = ImageDraw.Draw(image_summary)
            if platform.system() == "Darwin":
                pretty_font = ImageFont.truetype("Keyboard.ttf", 80)
            else:
                try:
                    pretty_font = ImageFont.truetype("arial.ttf", 80)
                except:
                    try:
                        pretty_font = ImageFont.load_default(size=80)
                    except:
                        pretty_font = ImageFont.load_default()
            title.text((20, 20),
                       "{}: {} analysis".format(name, table_name),
                       font=pretty_font,
                       fill=(0,0,0))
            image_summary.save("summary-{}-{}.png".format(name, table_name))
            if show:
                image_summary.show()
            # cleaning up
            if cleanup:
                for img_path in images_path:
                    os.remove(img_path)
                

    def ccz_identifier(self):
        result = ""
        for a in sorted(["LAT", "DDT", "DEG", "THI"]):
            result += pretty_spectrum(self.spectra[a])
        return result
