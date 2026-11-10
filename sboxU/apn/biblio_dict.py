# Maps a (n, source_index) pair to a bibliographic reference string.
# The source_index stored in an APNFunctions_compact entry is a key into
# the per-dimension dict below.

biblio_6 = {
    -1: "NOT SPECIFIED",
    0: "Banff",
    1: "Edel & Pott",
}

biblio_7 = {
    -1: "NOT SPECIFIED",
}

biblio_8 = {
    -1: "NOT SPECIFIED",
    0: "TO SOURCE NICE POLYNOMIALS LIST",
    1: "Yu, Y., Wang, M. & Li, Y. A matrix approach for constructing quadratic APN functions. Des. Codes Cryptogr. 73, 587–600 (2014). https://doi.org/10.1007/s10623-014-9955-3",
    2: "Yu, Y., Perrin, L. Constructing more quadratic APN functions with the QAM method. Cryptogr. Commun. 14, 1359–1369 (2022). https://doi.org/10.1007/s12095-022-00598-z",
    3: "C. Beierle and G. Leander, New Instances of Quadratic APN Functions, in IEEE Transactions on Information Theory, vol. 68, no. 1, pp. 670-678, Jan. 2022, doi: 10.1109/TIT.2021.3120698",
    4: "Beierle, C., Leander, G., & Perrin, L. (2022). Quadratic APN Extensions (2.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5821018",
    5: "Weng, G., Tan, Y., & Gong, G. (2013). On quadratic almost perfect nonlinear functions and their related algebraic object. In Workshop on Coding and Cryptography, WCC.",
    6: "Beierle, C., Langevin, P., Leander, G., Polujan, A., & Rasoolzadeh, S. (2025). Millions of inequivalent quadratic APN functions in eight variables [Data set]. Zenodo. https://doi.org/10.5281/zenodo.16752428",
}

biblio_9 = {
    -1: "NOT SPECIFIED",
    0: "Gold Monomial",
    1: "Kasami Monomial",
    2: "Welch Monomial",
    3: "Inverse Monomial",
    4: "L. Budaghyan, C. Carlet, G. Leander. Constructing new APN functions from known ones.",
    5: "L. Budaghyan, C. Carlet, G. Leander. On a construction of quadratic APN functions.",
    6: "L. Budaghyan, M. Calderini, C. Carlet, R. S. Coulter, I. Villa. Constructing APN Functions through Isotopic Shifts",
    7: "L. Budaghyan, C. Carlet, G. Leander. Constructing new APN functions from known ones",
    8: "LK23",
    9: "LZLQ22a"
}

biblio_10 = {
    -1: "NOT SPECIFIED",
    0: "Gold Monomial",
    1: "Kasami Monomial",
    2: "Dobbertin Monomial",
    3: "BCL09",
    4: "T19",
    5: "BHK20",
    6: "BCV20",
    7: "G22",
    8: "CLV22",
    9: "ZKLPT22",
    10: "LZLQ22b",
    11: "SPZ25",
    12: "TBD"
}

_BIBLIO = {6: biblio_6, 7: biblio_7, 8: biblio_8, 9: biblio_9, 10: biblio_10}


def from_which_paper(n, k):
    """Returns the bibliographic reference for source index `k` among n-bit APN functions,
    or "Not sourced yet" if no entry is registered for that (n, k) pair."""
    d = _BIBLIO.get(n, {})
    return d.get(k, "Not sourced yet")
