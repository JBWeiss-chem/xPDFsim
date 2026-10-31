"""X-ray and neutron pair distribution functions.

The definitions follow P. F. Peterson, D. Olds, M. T. McDonnell & K. Page,
J. Appl. Cryst. 54, 317-332 (2021) and the corrigendum by P. F. Peterson &
D. A. Keen, J. Appl. Cryst. 54, 1542-1545 (2021). For neutrons, b_i is the
coherent scattering length of the element (Neutron Data Booklet, Sears'
values). For X-rays, it is approximated by the X-ray scattering length at
Q = 0, b_i = Z_i r_e, i.e. the Q dependence of the atomic form factors is
neglected.

The functions are passed around as tables: dictionaries of numpy arrays with
the key "radii", one key per element pair (the partials) and "total".
"""

import itertools
import numpy as np

from .elements import ATOMIC_NUMBERS, NEUTRON_B

R_E = 2.8179403262      # classical electron radius in fm
FM2_PER_BARN = 100.0
RADIATIONS = ("xray", "neutron")

# Real-space functions of Table 4 of Peterson et al. (2021) with their units,
# the most commonly used first
CONVENTIONS = {
    "G_r": "Å⁻²",                # reduced pair distribution function, eq. (19)
    "g_r": "",                   # pair distribution function, eq. (22)
    "D_r": "barn Å⁻²",           # differential correlation function, eq. (28)
    "GK_r": "barn",              # total radial distribution function, eq. (26)
    "R_r": "Å⁻¹",                # weighted radial distribution function, eq. (15)
    "rho_r": "Å⁻³",              # density function, eq. (16)
    "T_r": "barn Å⁻²",           # total correlation function, eq. (29)
    "N_r": "",                   # coordination number, eq. (31)
}


def gen_R_r(structure, r_max, dr, verbose=True, radiation="xray", r_extra=1.0):
    """Takes a Structure, calculates neighbor lists, transforms them
    into a histogram, and finally transforms them into the weighted radial
    distribution function R(r) [Peterson et al. (2021), eq. (15)].

    Args:
        structure (Structure): Structure of which the PDF shall be calculated.
        r_max (float or int): Maximum r to which the PDF shall be calculated.
        dr (float or int): Bin width of the histogram.
        verbose (Boolean): If True, progress messages are printed.
        radiation (string): "xray" or "neutron", sets the scattering lengths.
        r_extra (float): The histograms go on to r_max + r_extra, so that
        broadening near r_max includes the peaks beyond it (about 6 sigma).

    Returns:
        dictionary: Table of the R(r) histograms. The partials are the
        weighted contributions of the element pairs, which add up to the total.
        float: Number density of the material.
    """
    # composition with partial occupancies counted fractionally
    elements = list(structure.composition())
    total_atoms = structure.num_atoms()
    rho = total_atoms / structure.volume
    weights, _ = pair_weights(structure, radiation)

    # Histogram bins; calculated up to r_max + r_extra so that broadening near
    # r_max is not truncated
    n_bins = int(round((r_max + r_extra) / dr))
    radii = (np.arange(n_bins) + 0.5) * dr

    if verbose:
        if structure.is_ordered:
            print("Ordered structure detected.\n")
        else:
            print("Disordered structure detected.\n")
        print("Calculating PDF...\n")

    occupancies = get_occupancies(structure, elements)
    centers, neighbors, _, distances = structure.neighbor_list(r_max + r_extra)

    # Histogram distances of every element pair, weighted by the occupancies
    # of both sites (1 for ordered structures); the bins are uniform, so the
    # bin of every distance is computed once and counted with bincount
    bins = np.floor(distances / dr).astype(np.int64)
    inside = bins < n_bins
    bins, centers, neighbors = bins[inside], centers[inside], neighbors[inside]
    table = {"radii": radii}
    for i, j in itertools.combinations_with_replacement(range(len(elements)), 2):
        occupancy_weights = occupancies[centers, i] * occupancies[neighbors, j]
        if i != j:
            occupancy_weights += occupancies[centers, j] * occupancies[neighbors, i]
        hist = np.bincount(bins, weights=occupancy_weights, minlength=n_bins)

        label = pair_label(elements[i], elements[j])
        table[label] = weights[label]["bb"] / (total_atoms * dr) * hist

    table["total"] = sum(table[label] for label in partial_labels(table))
    return table, rho


def partial_labels(table):
    """The element pairs of a table."""
    return [label for label in table if label not in ("radii", "total")]


def pair_label(element_1, element_2):
    """Label of an element pair, e.g. 'Cl-Na'."""
    return "-".join(sorted([element_1, element_2]))


def scattering_length(element, radiation):
    """Scattering length b of an element in fm.

    Args:
        element (string): Element symbol.
        radiation (string): "xray" (b = Z r_e) or "neutron" (coherent
        scattering length).

    Returns:
        float: Scattering length in fm.
    """
    if radiation == "xray":
        return ATOMIC_NUMBERS[element] * R_E
    if radiation == "neutron":
        try:
            return NEUTRON_B[element]
        except KeyError:
            raise ValueError(f"No neutron scattering length is tabulated for {element}.")
    raise ValueError(f"Unknown radiation {radiation}")


def pair_weights(structure, radiation="xray"):
    """Weights of all element pairs for X-rays or neutrons.

    Args:
        structure (Structure): Structure of which the PDF shall be calculated.
        radiation (string): "xray" or "neutron".

    Returns:
        dictionary: For every pair label, "bb" = b_i b_j / <b_coh>^2 (the factor
        in eq. 15) and "W" = the Faber-Ziman weight c_i c_j b_i b_j / <b_coh>^2
        of eq. (61), summed over both orders of unlike pairs, so that the
        weights of all pairs add up to 1. With neutrons, elements with a
        negative scattering length (e.g. H, Ti, Mn) give negative weights.
        float: <b_coh>^2 in barn (eq. 57).

    Raises:
        ValueError: If a scattering length is missing or <b_coh> is zero.
    """
    composition = structure.composition()
    total = sum(composition.values())
    fraction = {el: amount / total for el, amount in composition.items()}
    elements = list(composition)
    b = {el: scattering_length(el, radiation) for el in elements}
    mean_b = sum(fraction[el] * b[el] for el in elements)
    if abs(mean_b) < 1e-6:
        raise ValueError("The average scattering length is zero, so the normalised functions are not defined.")
    weights = {}
    for el_1, el_2 in itertools.combinations_with_replacement(elements, 2):
        bb = b[el_1] * b[el_2] / mean_b**2
        c_1, c_2 = fraction[el_1], fraction[el_2]
        multiplicity = 1 if el_1 == el_2 else 2
        weights[pair_label(el_1, el_2)] = {"bb": bb, "W": multiplicity * c_1 * c_2 * bb}
    mean_b2 = mean_b**2 / FM2_PER_BARN
    return weights, mean_b2


def linear_terms(convention, r, rho, mean_b2):
    """Every function X(r) of Table 4 except N(r) is linear in R(r):
    X(r) = a(r) R(r) + c(r) W, with W = 1 for the total and the Faber-Ziman
    weight for the partials (Appendix A1 and F of Peterson et al., 2021).

    Returns:
        numpy array: a(r).
        numpy array: c(r).
    """
    zero = np.zeros_like(r)
    if convention == "R_r":
        return np.ones_like(r), zero
    if convention == "rho_r":
        return 1 / (4 * np.pi * r**2), zero
    if convention == "g_r":
        return 1 / (4 * np.pi * rho * r**2), zero
    if convention == "GK_r":
        return mean_b2 / (4 * np.pi * rho * r**2), zero - mean_b2
    if convention == "G_r":
        return 1 / r, -4 * np.pi * rho * r
    if convention == "D_r":
        return mean_b2 / r, -4 * np.pi * rho * mean_b2 * r
    if convention == "T_r":
        return mean_b2 / r, zero
    raise ValueError(f"Unknown convention {convention}")


def convert_R_r(table, rho, convention, weights, mean_b2):
    """Converts R(r) histograms into another real-space function.

    Args:
        table (dictionary): Table generated by gen_R_r().
        rho (float): Number density of the material.
        convention (string): One of CONVENTIONS.
        weights (dictionary): Pair weights generated by pair_weights().
        mean_b2 (float): <b_coh>^2 in barn generated by pair_weights().

    Returns:
        dictionary: Table of the converted functions. N(r) is given as total
        only, as it is not weighted by the scattering lengths.
    """
    radii = table["radii"]
    labels = partial_labels(table)
    out = {"radii": radii}

    if convention == "N_r":
        # Unweighted RDF(r) of eq. (14), integrated up to the bin centres
        dr = radii[1] - radii[0]
        rdf = sum(table[label] / weights[label]["bb"] for label in labels)
        out["total"] = (np.cumsum(rdf) - rdf / 2) * dr
        return out

    a, c = linear_terms(convention, radii, rho, mean_b2)
    for label in labels:
        out[label] = a * table[label] + c * weights[label]["W"]
    out["total"] = a * table["total"] + c
    return out


def baselines(convention, r, rho, weights, mean_b2, labels):
    """The terms c(r) W of linear_terms() on the grid r, for every label of a table.

    They are smooth (constant or linear in r), so broadening must leave them
    unchanged; broadening them with the histograms would bend them where the
    Gaussians are cut off, at r = 0 and at the end of the histograms.

    Returns:
        dictionary: Label -> c(r) W (zero for N(r)).
    """
    if convention == "N_r":
        return {label: np.zeros_like(r) for label in labels}
    with np.errstate(divide="ignore"):         # a(r) is infinite at r = 0, c(r) is not
        _, c = linear_terms(convention, r, rho, mean_b2)
    return {label: c *(1.0 if label == "total" else weights[label]["W"]) for label in labels}


def write_csv(table, filename):
    """Writes a table as a CSV file: radii, partials, total."""
    labels = ["radii"] + [label for label in table if label != "radii"]
    data = np.column_stack([table[label] for label in labels])
    np.savetxt(filename, data, delimiter=",", header=",".join(labels), comments="", fmt="%.10g")


def asymptotes(convention, r, rho, mean_b2):
    """Low-r and high-r limits of Table 4 of Peterson et al. (2021).

    Returns:
        list: The limits as arrays over r.
    """
    zero = np.zeros_like(r)
    return {
        "R_r": [zero, 4 * np.pi * rho * r**2],
        "rho_r": [zero, zero + rho],
        "g_r": [zero, zero + 1],
        "GK_r": [zero - mean_b2, zero],
        "G_r": [-4 * np.pi * rho * r, zero],
        "D_r": [-4 * np.pi * rho * mean_b2 * r, zero],
        "T_r": [zero, 4 * np.pi * rho * mean_b2 * r],
        "N_r": [zero, 4 / 3 * np.pi * rho * r**3],
    }[convention]


def get_occupancies(structure, elements):
    """Builds a matrix containing the occupancy of every element on every site.

    Args:
        structure (Structure): Structure of which the PDF shall be calculated.
        elements (list): Element symbols present in the structure.

    Returns:
        numpy array: Occupancies with shape (number of sites, number of
        elements).
    """
    element_index = {el: k for k, el in enumerate(elements)}
    occupancies = np.zeros((len(structure), len(elements)))
    for site_index, species in enumerate(structure.species):
        for el, occupancy in species.items():
            occupancies[site_index, element_index[el]] += occupancy
    return occupancies
