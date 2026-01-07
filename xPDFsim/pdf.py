import numpy as np
import pandas as pd
from pymatgen.core import Element
import itertools as iter
from .utils import *


def gen_R_r(structure, r_max, dr):
    """Takes a pymatgen Structure, calculates neighbor lists, transforms them 
    into a histogram, and finally transforms them into X-ray pair distribution 
    functions.

    Args:
        structure (pymatgen Stcucture): Structure of which the PDF shall be 
        calculated.
        r_max (float or int): Maximum r to which the PDF shall be calculated.
        dr (float or int): Bin width of the histogram.

    Returns:
        dataframe: Histogramms of the R(r) functions.
        float: Number density of the material.
    """
    # Setup for later
    total_atoms, c_is = get_fractions(structure.formula)
    rho = total_atoms / structure.lattice.volume

    # Generate element pairs for partial PDFs
    element_list = [element.symbol for element in structure.elements]
    element_pairs_raw = list(iter.combinations_with_replacement(element_list, 2))
    sorted_elements = [sorted(tup) for tup in element_pairs_raw]
    element_pairs = [f'{el[0]}-{el[1]}' for el in sorted_elements]
    
    # Calculate X-ray normalisation factor
    K_is = [Element(atom_type).Z for atom_type in element_list]
    norm_x_ray = 0
    for i in range(len(c_is)):
        norm_x_ray += c_is[i] * K_is[i]
    norm_x_ray = norm_x_ray**2 

    radii = np.arange(dr, r_max+1, dr) 

    # Check if structure is ordered to safe some computing time
    if structure.is_ordered:
        print("Ordered structure detected.\n")
        print("Calculating PDF...\n")
        rdf_sharp = gen_sharp_rdf(structure, element_pairs, r_max)
        hist_dict = gen_histogram(rdf_sharp, radii)

    elif structure.is_ordered == False:
        print("Disordered structure detected.\n")
        print("Calculating PDF...\n")
        rdf_sharp, rdf_sharp_weights = gen_sharp_rdf_dis(structure, 
                                                         element_pairs, r_max)
        hist_dict = gen_histogram_dis(rdf_sharp, rdf_sharp_weights, radii)

    hist_dict['radii'] = radii[:-1]
    
    # Calculate R(r) from distance histograms
    R_r_df = gen_R_r_df(hist_dict, norm_x_ray, total_atoms, dr)
    
    return R_r_df, rho
    

def gen_sharp_rdf(structure, element_pairs, r_max):
    """Generates neighbor list and sorts distances by element pairs.

    Args:
        structure (pymatgen Stcucture): Structure of which the PDF shall be 
        calculated.
        element_pairs (list): Element-element pairs of which the PDF shall be 
        calculated.
        r_max (float or int): Maximum r to which the PDF shall be calculated.

    Returns:
        dictionary: Contains all Element-Element PDFs.
    """
    all_neighbors = structure.get_all_neighbors(r=r_max+1)
    rdf_sharp = {}
    for element_pair in element_pairs:
        rdf_sharp[element_pair] = []
    for id, site in enumerate(structure.sites):
        element_1 = str(site.specie)
        neighbors = all_neighbors[id]
        for neighbor in neighbors:
            element_2 = str(neighbor.specie)
            pair_list = sorted([element_1, element_2])
            pair_label = f'{pair_list[0]}-{pair_list[1]}'
            dist = np.linalg.norm(site.coords - neighbor.coords)
            rdf_sharp[pair_label].append(dist)
    return rdf_sharp


def gen_sharp_rdf_dis(structure, element_pairs, r_max):
    """Generates neighbor list and sorts distances by element pairs. For d
    isordered materials.

    Args:
        structure (pymatgen Stcucture): Structure of which the PDF shall be 
        calculated.
        element_pairs (list): Element-element pairs of which the PDF shall be 
        calculated.
        r_max (float or int): Maximum r to which the PDF shall be calculated.

    Returns:
        dictionary: Contains all Element-Element PDFs.
    """    
    all_neighbors = structure.get_all_neighbors(r=r_max+1)
    rdf_sharp = {}
    rdf_sharp_weights = {}
    for element_pair in element_pairs:
        rdf_sharp[element_pair] = []
        rdf_sharp_weights[element_pair] = []
    for id, site in enumerate(structure.sites):
        species_1 = site.species_string.split(':')
        element_1 = species_1[0]
        if len(species_1) == 1:
            occu_1 = 1
        elif len(species_1) == 2:
            occu_1 = float(species_1[1])
        neighbors = all_neighbors[id]
        for neighbor in neighbors:
            species_2 = neighbor.species_string.split(':')
            element_2 = species_2[0]
            if len(species_2) == 1:
                occu_2 = 1
            elif len(species_2) == 2:
                occu_2 = float(species_2[1])
            pair_list = sorted([element_1, element_2])
            pair_label = f'{pair_list[0]}-{pair_list[1]}'
            dist = np.linalg.norm(site.coords - neighbor.coords)
            weight = occu_1 * occu_2
            rdf_sharp[pair_label].append(dist)
            rdf_sharp_weights[pair_label].append(weight)
    return rdf_sharp, rdf_sharp_weights


def gen_histogram(rdf_dict, radii):
    """Turns distance lists into histograms.

    Args:
        rdf_dict (dictionary): Contains lists of exact distances.
        radii (numpy array): Radii corresponding to the to-be-generated 
        histogram.

    Returns:
        dictionary: Contains arrays of the generated histograms.
    """
    hist_dict = {}
    for pair in rdf_dict:
        hist_dict[pair] = np.histogram(rdf_dict[pair], bins=radii)[0]
    return hist_dict


def gen_histogram_dis(rdf_dict, rdf_dict_weights, radii):
    """Turns distance lists into histograms. For disordered materials.

    Args:
        rdf_dict (dictionary): Contains lists of exact distances.
        rdf_dict_weights (dictionary): Contains weights for each distance 
        calculated from occupancies.
        radii (numpy array): Radii corresponding to the to-be-generated 
        histogram.

    Returns:
        dictionary: Contains arrays of the generated histograms.
    """

    hist_dict = {}
    for pair in rdf_dict:
        hist_dict[pair] = np.histogram(rdf_dict[pair], bins=radii, 
                                       weights=rdf_dict_weights[pair])[0]
    return hist_dict


def gen_R_r_df(hist_dict, norm_x_ray, total_atoms, dr):
    """Adds X-ray contributions to histograms of distances.

    Args:
        hist_dict (dictionary): Contains histograms of distances.
        norm_x_ray (float): Square of average atomic number in the structure.
        total_atoms (int): Total atoms in the unit cell.
        dr (float or int): Bin width of the histogram.

    Returns:
        dataframe: Contains R(r) functions.
    """
    R_r_dict = {}
    for pair in hist_dict:
        if pair == 'radii':
            break
        element_1 = pair.split('-')[0]
        element_2 = pair.split('-')[1]
        K_1 = Element(element_1).Z
        K_2 = Element(element_2).Z
        R_r = K_1*K_2 / (norm_x_ray*total_atoms*dr) * hist_dict[pair]
        R_r_dict[pair] = R_r

    R_r_dict['radii'] = hist_dict['radii']
    df = pd.DataFrame(R_r_dict).round(4)
    df['total'] = df.drop(columns=['radii']).sum(axis=1)
    return df
