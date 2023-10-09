#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont

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
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
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
            result.append("         1:{}".format(pretty_vector(self.struct[c][1])))
        return result
    

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"

            
class CycleAnomaly:
    def __init__(self, s):
        self.name = "cycle decomposition"
        self.cycles = cycle_decomposition(s)
        self.spectrum = defaultdict(int)
        for c in self.cycles:
            self.spectrum[len(c)] += 1
        # !TODO! evaluate the anomaly of all cycle types 
        self.positive_anomaly = 0
        self.negative_anomaly = 0

    def summary(self):
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
    def __init__(self, s, threshold=40):
        self.name = "CCZ"
        N = int(log(len(s), 2))
        self.spaces = []
        self.spaces_all_found = True
        # doing the thickness spectrum by hand so that this does not
        # take up too much time a priori
        self.lat_zeroes = lat_zeroes(s)
        for space in vector_spaces_bases_iterator(self.lat_zeroes, N, 2*N):
            self.spaces.append(space)
            if len(self.spaces) == threshold:
                self.spaces_all_found = False
                break
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
    for a in range(2, len(s)):
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
                 # how to perform the analysis
                 store_tables=None,
                 deep=False,
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
        self.advanced = {}
        if deep:
            for table_name in sorted(self.tables.keys()):
                t = self.tables[table_name]
                threshold = 4 if table_name == "DDT" else 2
                row_spectra = sort_table_rows(t, threshold=threshold)
                col_spectra = sort_table_rows([[t[a][b] for a in range(0, 2**self.N)]
                                               for b in range(0, 2**self.N)],
                                              threshold=threshold)
                if len(row_spectra.keys()) > 0:
                    self.advanced["{} rows common spectra".format(table_name)] = row_spectra
                if len(col_spectra.keys()) > 0:
                    self.advanced["{} cols common spectra".format(table_name)] = col_spectra
                    
            # more sophisticated analysis
            # !TODO! advanced analysis

            
            
    def expected_distribution(self, table_name):
        if table_name == "DDT":
            return ddt_coeff_probability
        elif table_name == "LAT":
            if is_permutation(self.lut):
                return lat_coeff_probability_permutation
            else:
                return lat_coeff_probability_function
        elif table_name == "BCT":
                return bct_coeff_probability
        else:
            raise Exception("unknown table name: " + table_name)


    def show(self,
             indent="",    
             noteworthy_threshold=10,
             noteworthy_mark="  # Noteworty"):
        print(indent + "LUT:")
        print(indent + "      {}".format(
            pretty_vector(list(range(0, 16)), template="{:02x} ")[1:-1]
        ))
        print(indent + "-"*54)
        for i in range(0, 2**(self.N-4)):
            print(indent + "{:2x}0 | {}".format(
                i, pretty_vector(self.lut[i*16:i*16+16],
                                 template="{:02x} "
                                 )[1:-1]
            ))
        good_properties = []
        bad_properties  = []
        for anomaly_name in sorted(self.anomalies.keys()):
            a = self.anomalies[anomaly_name]        
            to_print = indent + a.name
            if a.positive_anomaly >= noteworthy_threshold:
                to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                good_properties.append(a.name)
            elif a.negative_anomaly >= noteworthy_threshold:
                to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                bad_properties.append(a.name)
            print("\n# {}\n".format(to_print))
            for line in a.summary():
                print(indent + line)
        if len(good_properties) > 0:
            print("\n" + indent + "# Noteworthy qualities: " + str(good_properties)[1:-1])
        if len(bad_properties) > 0:
            print("\n" + indent + "# Noteworthy flaws    : " + str(bad_properties)[1:-1])
        if (len(good_properties) == 0) and (len(bad_properties) == 0):
            print(indent + "The function is dull")
        if len(self.advanced.keys()) > 0:
            print(indent + "Advanced properties:")
            for name in self.advanced.keys():
                print(indent + " - " + name)
                print(indent + "   " + str(self.advanced[name]))


    def save_pollock(self, 
                     name=None,
                     vmins=[None, None], 
                     vmaxs=[None, None], 
                     color_schemes=["CMRmap_r", "seismic"],
                     show=True):
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
            images = []
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
                    color_scheme=color_schemes[i],
                )
                images.append(Image.open(file_name+".png"))
            # comparison with expected distribution
            expected_distrib = None
            if table_name == "DDT":
                expected_distrib = ddt_coeff_probability
            elif table_name == "LAT":
                if is_permutation(self.lut):
                    expected_distrib = lat_coeff_probability_permutation
                else:
                    expected_distrib = lat_coeff_probability_function
            elif table_name == "BCT":
                expected_distrib = bct_coeff_probability
            file_name = "distrib-{}-{}".format(name, table_name)
            plot_statistical(
                spec,
                n=self.N,
                expected_distrib=expected_distrib,
                file_name=file_name
            )
            images.append(Image.open(file_name+".png"))
            # patterns in the variance
            file_name = "{}-var-{}".format(name, table_name)
            plot_table_variances(self.tables[table_name], file_name=file_name)
            images.append(Image.open(file_name+".png"))
            # patterns by row and columns
            file_name = "distrib-by-row-{}-{}".format(name, table_name)
            plot_statistical_by_rows(
                self.tables[table_name],
                n=self.N,
                expected_distrib=self.expected_distribution(table_name),
                file_name=file_name
            )
            images.append(Image.open(file_name+".png"))
            file_name = "distrib-by-col-{}-{}".format(name, table_name)
            plot_statistical_by_rows(
                [[self.tables[table_name][a][b] for a in range(0, 2**self.N)]
                 for b in range(0, 2**self.N)],
                n=self.N,
                expected_distrib=self.expected_distribution(table_name),
                file_name=file_name
            )
            images.append(Image.open(file_name+".png"))
            # assembling summary picture
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
            title.text((20, 20),
                       "{}: {} analysis".format(name, table_name),
                       font=ImageFont.truetype("arial.ttf", 80),
                       fill=(0,0,0))
            image_summary.save("summary-{}-{}.png".format(name, table_name))
            if show:
                image_summary.show()


    def ccz_identifier(self):
        result = ""
        for a in sorted(["LAT", "DDT", "DEG", "THI"]):
            result += pretty_spectrum(self.spectra[a])
        return result
