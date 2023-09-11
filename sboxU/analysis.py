#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint

from collections import defaultdict

from .sboxU_cython import *
from .utils import *
from .display import *
from .diff_lin import *
from .ccz import *


# !SECTION! Statistical anomalies 
# -------------------------------

class differential_anomaly:
    """The differential anomaly quantifies the probability that a
    random function has a pair (u, n(u)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, u is the differential uniformity, and n(u) is
    its number of occurrences in the DDT.

    """
    def __init__(self, s):
        self.name = "differential"
        self.spectrum = differential_spectrum(s)
        self.positive_anomaly = table_anomaly(s, "DDT", spec=self.spectrum)
        self.negative_anomaly = table_negative_anomaly(s, "DDT", spec=self.spectrum)

    def summary(self):
        return [
            "spectrum:   {}".format(pretty_spectrum(self.spectrum)),
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
        ]

    def __str__(self):
        for line in self.summary():
            result += line + "\n"
        
    
class linear_anomaly:
    """The linear anomaly quantifies the probability that a
    random function has a pair (l, n(l)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, l is the linearity, and n(l) is
    the sum of the number of occurrences of l and -l in the LAT.

    """
    def __init__(self, s):
        self.name = "linear"
        self.spectrum = walsh_spectrum(s)
        self.positive_anomaly = table_anomaly(s, "LAT", spec=self.spectrum)
        self.negative_anomaly = table_negative_anomaly(s, "LAT", spec=self.spectrum)

    def summary(self):
        return [
            "spectrum:   {}".format(pretty_spectrum(self.spectrum, absolute=True)),
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
        ]

    def __str__(self):
        for line in self.summary():
            result += line + "\n"

    
class boomerang_anomaly:
    """The boomerang anomaly quantifies the probability that a
    random function has a pair (b, n(b)) that is lexicographically
    lower (positive anomaly) or higher (negative anomaly) than that of
    the studied function.

    In this expression, b is the boomerang uniformity, and n(b) is
    its number of occurrences in the BCT.

    """
    def __init__(self, s):
        self.name = "boomerang"
        self.spectrum = boomerang_spectrum(s)
        self.positive_anomaly = table_anomaly(s, "BCT", spec=self.spectrum)
        self.negative_anomaly = table_negative_anomaly(s, "BCT", spec=self.spectrum)

    def summary(self):
        return [
            "spectrum:   {}".format(pretty_spectrum(self.spectrum, absolute=True)),
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
        ]

    def __str__(self):
        for line in self.summary():
            result += line + "\n"


    
# !SECTION! Structural anomalies
# ------------------------------

class ccz_anomaly:
    def __init__(self, s):
        self.name = "CCZ"
        N = int(log(len(s), 2))
        self.spaces = get_lat_zeroes_spaces(s)
        self.spectrum = thickness_spectrum(s, spaces=self.spaces)
        self.negative_anomaly = 0
        self.positive_anomaly = 0
        for t in self.spectrum.keys():
            if t > 0:
                log_tu_density = N*2**N - 2**(N-t) * sum(lg2(i) for i in range(1, 2**t+1)) - t * (N-t) * 2**(N-t)
                log_affine_card = N + sum(lg2(2**N - 2**i) for i in range(0, N))
                self.positive_anomaly += (log_tu_density - 2*log_affine_card - N**2) * self.spectrum[t]

    def summary(self):
        return [
            "|Walsh zero spaces| = {}".format(len(self.spaces)),
            "thickness spectrum = {}".format(pretty_spectrum(self.spectrum)),
            "positive A: {:10.4f}".format(self.positive_anomaly),
            "negative A: {:10.4f}".format(self.negative_anomaly),
        ]

    def __str__(self):
        result = ""
        for line in self.summary():
            result += line + "\n"

    
class commutative_anomaly:
    def __init__(self, s):
        self.name = "commutative"
        N = int(log(len(s), 2))
        self.commutants = self_affine_equivalent_mappings(s)

        

# !SECTION! Automated analysis 
# ----------------------------

class Analysis:
    def __init__(self,
                 s,
                 name="sbox",
                 differential=True,
                 linear=True,
                 boomerang=True,
                 ccz=True,
                 commutative=False,
                 store_tables=True,
                 noteworthy_threshold=10,
                 noteworthy_mark="  # Noteworty",
                 indent="",
                 verbose=True):
        self.N = int(lg2(len(s)))
        self.name = name
        self.anomalies = {}
        self.tables = {}
        self.spectra = {}
        if differential:
            self.anomalies["differential"] = differential_anomaly(s)
            self.spectra["differential"] = self.anomalies["differential"].spectrum
            if store_tables:
                self.tables["DDT"] = ddt(s)
        if linear:
            self.anomalies["linear"] = linear_anomaly(s)
            self.spectra["linear"] = self.anomalies["linear"].spectrum
            if store_tables:
                self.tables["LAT"] = lat(s)
        if boomerang:
            self.anomalies["boomerang"] = boomerang_anomaly(s)
            self.spectra["boomerang"] = self.anomalies["boomerang"].spectrum
            if store_tables:
                self.tables["BCT"] = bct(s)
        if ccz:
            self.anomalies["CCZ"] = ccz_anomaly(s)
            self.spectra["thickness"] = self.anomalies["CCZ"].spectrum
        self.print()


    def print(self):
        good_properties = []
        bad_properties  = []
        if verbose:
            for anomaly_name in sorted(self.anomalies.keys()):
                a = self.anomalies[anomaly_name]        
                to_print = indent + a.name
                if a.positive_anomaly >= noteworthy_threshold:
                    to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                    good_properties.append(a.name)
                elif a.negative_anomaly >= noteworthy_threshold:
                    to_print = to_print + " "*(60-len(to_print)-len(noteworthy_mark)) + noteworthy_mark
                    bad_properties.append(a.name)
                print(to_print)
                for line in a.summary():
                    print(indent + " - " + line)
            if len(good_properties) > 0:
                print("Noteworthy qualities: " + str(good_properties)[1:-1])
            if len(bad_properties) > 0:
                print("Noteworthy flaws    : " + str(bad_properties)[1:-1])
            if (len(good_properties) == 0) and (len(bad_properties) == 0):
                print("The function is dull")


    def save_pollock(self, name=None):
        if name == None:
            name = self.name
        for table_name in self.tables.keys():
            save_pollock(
                self.tables[table_name], 
                name="{}-{}".format(name, table_name)
            )
