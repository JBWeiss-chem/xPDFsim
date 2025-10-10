import os
import numpy as np
from scipy.spatial import KDTree
from bokeh.plotting import figure, show
from bokeh.palettes import Category10_10
import argparse
import sys
import pandas as pd
from pymatgen.core import Structure, Element
import re


def gen_partial_rdf(sites_element_prim, sites_element_sc, shift_vec, radii, dr, eps):
    """
    Returns radial distribution function between two atom types.
    Atom type 1 is given as a primitive cell and atom type 2 as a supercell surrounding the primitive cell.
    """
    rdf = np.zeros(shape=(len(radii)))
    for occu_sc in sites_element_sc:
        
        coords_sc = [site["xyz"] for site in sites_element_sc[occu_sc]]
        
        # Create a KDTree for fast nearest-neighbor lookup of particles
        tree = KDTree(coords_sc)
        
        for r_idx, r in enumerate(radii):
            # compute n_i(r) for each particle in partial unit cell
            for occu_prim in sites_element_prim:
                mean_occu = (occu_sc+occu_prim) / 2
                coords_prim = [site["xyz"] + shift_vec for site in sites_element_prim[occu_prim]]
                for particle in coords_prim:
                    n = tree.query_ball_point(particle, r+dr-eps, return_length=True) - tree.query_ball_point(particle, r, return_length=True)
                    rdf[r_idx] += n*mean_occu
    return rdf


def get_fractions(formula):
    """
    Counts the total number of elements and calculates the fraction of each element
    in a chemical formula given as a string.
    """
    # Remove any whitespace from the formula
    formula = formula.replace(' ', '')

    # Regex to find element symbols and their counts
    pattern = r'([A-Z][a-z]*)(\d*\.?\d*)'
    elements = re.findall(pattern, formula)

    element_counts = {}
    total_atoms = 0.0

    for element, count in elements:
        # If count is not specified, it's 1
        if count == '':
            count = 1
        else:
            count = float(count)

        element_counts[element] = element_counts.get(element, 0) + count
        total_atoms += count

    element_fractions = {}
    for element, count in element_counts.items():
        element_fractions[element] = count / total_atoms
    
    fractions_list = list(element_fractions.values())

    return total_atoms, fractions_list


def find_equivalent_pairs(bond_list):
    """
    Finds equivalent pairs in a list of atom-atom corellation pairs (e.g. Cu–N and N-Cu).
    """
    # Dictionary to store the an atom pair as the key
    # and a list of its indices as the value.
    # Example: {'(N, Cu)': [2, 8]}
    pair_map = {}

    # Iterate through the input list to get both the index and the value
    for index, bond_str in enumerate(bond_list):
        # Split the string by the hyphen to get the two elements
        elements = bond_str.split('-')

        # Sort the elements alphabetically to create a canonical form.
        elements.sort()

        # Convert the sorted list to a tuple so it can be used as a dictionary key.
        canonical_pair = tuple(elements)

        if canonical_pair in pair_map:
            pair_map[canonical_pair].append(index)
        else:
            pair_map[canonical_pair] = [index]

    # Filter the results to include only the groups that have duplicates.
    result = [indices for indices in pair_map.values() if len(indices) > 1]

    return result
    
    
def find_intra_element_pairs(bond_list):
    """
    Finds indices of bonds between two of the same chemical element.
    """
    same_element_indices = []
    # Iterate through the list with an index
    for index, bond_str in enumerate(bond_list):
        # Split the bond string into its constituent elements
        elements = bond_str.split('-')
        # Check if the two elements are identical
        if elements[0] == elements[1]:
            # If they are, add the index to our list
            same_element_indices.append(index)
            
    return same_element_indices
    
    
def clean_g_rs(g_rs, g_r_labels):
    """
    Removes duplicate partial rdfs from a list of partial rdfs (e.g. Cu–N and N-Cu).
    """
    # Start by adding all intra-element g_rs to the list
    g_rs_clean = []
    g_r_labels_clean = []
    
    intra_element_pairs = find_intra_element_pairs(g_r_labels)
    for intra_element_pair in intra_element_pairs:
        g_rs_clean.append(g_rs[intra_element_pair])
        g_r_labels_clean.append(g_r_labels[intra_element_pair])

    # Find all duplicate pairs, add them together, and add them to new clean list
    equivalent_pairs = find_equivalent_pairs(g_r_labels)
    for equivalent_pair in equivalent_pairs:
        g_r_added = g_rs[equivalent_pair[0]] + g_rs[equivalent_pair[1]]
        g_rs_clean.append(g_r_added)
        g_r_labels_clean.append(g_r_labels[equivalent_pair[0]])

    return np.array(g_rs_clean), g_r_labels_clean


def build_sc(r_max, lattice):
    """
    Builds supercell with an appropriate size (according to r_max) 
    around unit cell. 
    """
    # get cell vector component of each direction
    x = lattice[0][0]
    y = lattice[1][1]
    z = lattice[2][2]
    
    # calculate how much one needs to extend in each direction
    sc_x = np.floor_divide(r_max, x*0.8) + 1
    sc_y = np.floor_divide(r_max, y*0.8) + 1
    sc_z = np.floor_divide(r_max, z*0.8) + 1

    # Create supercell matrix
    M = [[2*sc_x+1, 0, 0], 
         [0, 2*sc_y+1, 0], 
         [0, 0, 2*sc_z+1]]

    # calculate the vector for shifting the unit cell into the center of the supercell
    shift_vec = np.zeros(3)
    shift_vec += lattice[0] * sc_x
    shift_vec += lattice[1] * sc_y
    shift_vec += lattice[2] * sc_z

    return M, shift_vec


def broaden_pdfs(df, sigma):
    """
    Turns rdf histograms into broadened functions
    """
    # Create an array of x-values for the broadened spectrum
    radii = df["radii"].to_numpy()
    dr = radii[1] -radii[0]
    labels = df.columns.to_numpy()[1:]
    r_max = np.round(radii[-1]).astype(int) - 1
    radii_broad = np.linspace(0, r_max, r_max*50)
    norm_factor = dr / (sigma * np.sqrt(2 * np.pi))

    # Put everything into a dataframe
    data = {'radii': radii_broad}
    df_out = pd.DataFrame(data)

    for i, label in enumerate(labels):
        pdf = df.iloc[:, i+1].to_numpy()
        pdf_broad = np.zeros_like(radii_broad)
        # Apply Gaussian broadening for each mode
        for r, y in zip(radii, pdf):
            pdf_broad += y * norm_factor * np.exp(-((radii_broad - r)**2) / (2 * sigma**2))
        df_out[label] = pdf_broad
    return df_out.round(4)


def get_R_r(structure, r_max, dr, eps=1e-15, sigma=0.2):
    """
    This code builds a supercell around the unit cell of the provided pymatgen structure object.
    Element-element radial distribution functions are generated by iterating over
    the atoms in the unit cell and finding the corellarions with all atoms surrounding 
    them in a certain radius. The pair distribution functions are generated from the 
    radial distribution functions by taking the atomic number of each element into 
    account.

    Parameters
    ----------
    atoms : ase Atoms object
        The xPDF of this object shall be simulated
    n_sc : Integer
         Used to build a n_sc x n_sc x n_sc supercell
         this shall be swapped by r_max for the pdf generation
    dr : float
        Delta r. Determines the spacing between successive radii over which g(r)
        is computed.
    eps : float, optional
        Epsilon value used to find particles less than or equal to a distance 
        in KDTree.
    """

    sdict = structure.as_dict()
    lattice = np.array(sdict["lattice"]["matrix"])
    sites = sdict["sites"]
    atom_types = [element.symbol for element in structure.elements]
    total_atoms, c_is = get_fractions(structure.formula)
    rho = total_atoms / sdict["lattice"]["volume"]


    # Calculate K_i values (atomic numbers) for all elements
    # correct this by average occupancy
    K_is = [Element(atom_type).Z for i, atom_type in enumerate(atom_types)]

    # calculate average atomic number for normalisation
    norm_x_ray = 0
    for i in range(len(c_is)):
        norm_x_ray += c_is[i] * K_is[i]
    norm_x_ray = norm_x_ray**2  

    # build supercell based on r_max
    M, shift_vec = build_sc(r_max+1, lattice)

    # Create a dictionary for each element sorted by occupancies
    sites_sorted_prim = []
    for atom_type in atom_types:
        element_sites = {}
        
        # Create list of occuring occupancies   
        occus = []
        for site in sites:
            if site["species"][0]["element"] == atom_type:
                if site["species"][0]["occu"] not in occus:
                    occus.append(site["species"][0]["occu"])

        for occu in occus:
            element_sites[occu] = []

        # Sort sites according to occupancy
        for site in sites:
            if site["species"][0]["element"] == atom_type:
                occu = site["species"][0]["occu"]
                element_sites[occu].append(site)
        sites_sorted_prim.append(element_sites)
    

    # For each occupancy of each element, build a supercell
    sites_sorted_sc = []
    for element_sites in sites_sorted_prim:
        sdict_sc = structure.as_dict()
        sites_element_sc = {}
        for occu in element_sites:
            sites_occu = element_sites[occu]
            sdict_sc["sites"] = sites_occu
            structure_sc = Structure.from_dict(sdict_sc)
            structure_sc.make_supercell(M)
            sdict_sc = structure_sc.as_dict()
            sites_sc = sdict_sc["sites"]
            sites_element_sc[occu] = sites_sc
        sites_sorted_sc.append(sites_element_sc)


    # setup for final output
    R_rs = []
    R_r_labels = []
    radii = np.arange(dr, r_max+1, dr) # radii for generation of all g_rs

    # Main function
    # iterate over primitive cells and ocalculate correlations with supercells
    for element_id_prim, sites_element_prim in enumerate(sites_sorted_prim):
        # set values for X-ray contributions to G_r
        K_i = K_is[element_id_prim]
        
        # now iterate over supercells for calculating correlations with primitive cell
        for element_id_sc, sites_element_sc in enumerate(sites_sorted_sc):
            
            # calculate rdf for this element pair
            rdf = gen_partial_rdf(sites_element_prim, sites_element_sc, shift_vec, radii, dr, eps)

            # Add X-ray contributions 
            K_j = K_is[element_id_sc]
            R_r = K_i*K_j / (norm_x_ray*total_atoms*dr) * rdf

            # add everything to lists of all partials
            R_rs.append(R_r)
            R_r_label = f"{atom_types[element_id_prim]}-{atom_types[element_id_sc]}"
            R_r_labels.append(R_r_label)

    # Remove duplicates
    R_rs_clean, R_r_labels_clean = clean_g_rs(R_rs, R_r_labels)

    # Put everything into a dataframe
    data = {'radii': radii}
    total_pdf = np.zeros_like(radii)
    for i, col_label in enumerate(R_r_labels_clean):
        data[col_label] = R_rs_clean[i]
        total_pdf += R_rs_clean[i]
    data["total"] = total_pdf
    df = pd.DataFrame(data).round(4)
    
    return df, rho
    

def format_plot(p):
    """
    Formats the Brokeh plot nicely
    """
    p.xaxis.axis_label_text_font_size = "20pt"
    p.yaxis.axis_label_text_font_size = "20pt"
    p.yaxis.major_label_text_font_size = "20pt"
    p.xaxis.major_label_text_font_size = "20pt"
    p.legend.label_text_font_size = "20pt"
    p.legend.glyph_width = 70

    # Enable interactive legend
    p.legend.click_policy="hide"


def plot_R_r(df, rho, sigma, histogram):
    """
    Turns dataframe into R(r) plot
    """
    # Apply gaussian broadening
    df_broad = broaden_pdfs(df, sigma)
    if histogram == True:
        df_broad = df
    radii = df_broad["radii"].to_numpy()
    labels = df.columns.to_numpy()[1:-1]

    R_r_total = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, x_axis_label='r / Å', y_axis_label="R(r)", sizing_mode="stretch_both")
    
    # Plot full PDF
    p.line(radii, R_r_total, color="black", line_width=1.5, legend_label="R(r)")

    for label, color in zip(labels, Category10_10):
        R_r = df_broad[label].to_numpy()
        p.line(radii, R_r, color=color,line_width=1.5, legend_label=label)

    format_plot(p)
    p.legend.location = "top_left"

    show(p)
    return df_broad


def plot_g_r(df, rho, sigma, histogram):
    """
    Turns dataframe into g(r) plot
    """
    radii = df["radii"].to_numpy()
    labels = df.columns.to_numpy()[1:]

    # Calculate g_r partials from R_r partials
    for label in labels:
        df[label] = df[label] / (4*np.pi*rho*radii**2)

    # Apply gaussian broadening
    df_broad = broaden_pdfs(df, sigma)
    if histogram == True:
        df_broad = df
    radii_broad = df_broad["radii"].to_numpy()
    g_r_total = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, x_axis_label='r / Å', y_axis_label="g(r)", sizing_mode="stretch_both")

    # Plot full PDF
    p.line(radii_broad, g_r_total, color="black", line_width=1.5, legend_label="g(r)")    

    for label, color in zip(labels[:-1], Category10_10):
        g_r = df_broad[label].to_numpy()
        p.line(radii_broad, g_r, color=color,line_width=1.5, legend_label=label)

    p.legend.location = "top_right"
    format_plot(p)
        
    show(p)
    return df_broad
    

def plot_G_r(df, rho, sigma, histogram):
    """
    Turns dataframe into G(r) plot
    """
    radii = df["radii"].to_numpy()
    total_R_r = df["total"].to_numpy()

    G_r = total_R_r / radii - 4 * np.pi * rho * radii
    data = {'radii': radii}
    data["total"] = G_r
    df_G_r = pd.DataFrame(data)
    
    # Apply gaussian broadening
    df_broad = broaden_pdfs(df_G_r, sigma)
    if histogram == True:
        df_broad = df_G_r
    radii_broad = df_broad["radii"].to_numpy()
    G_r_broad = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, x_axis_label='r / Å', y_axis_label="G(r)", sizing_mode="stretch_both")

    # Plot full PDF
    p.line(radii_broad, G_r_broad, color="black", line_width=1.5, legend_label="G(r)")

    p.legend.location = "top_right"
    format_plot(p)
    show(p)
    return df_broad


def main(args):
    print("")
    try:
        print(f"Input file received: {args.input_file}\n")
        structure = Structure.from_file(args.input_file, primitive=False)
        print("\nSuccessfully processed input file.\n")
    except FileNotFoundError:
        print(f"Error: The file '{args.input_file}' was not found.")
        sys.exit(1) # Exit with an error code
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
    
    print("Calculating PDF...\n")
    df, rho = get_R_r(structure, args.r_max, dr=args.bin_width)

    print("Plotting...\n")
    if args.pdf_type == "R_r":
        df_out = plot_R_r(df, rho, args.sigma, args.histogram)

    if args.pdf_type == "g_r":
        df_out = plot_g_r(df, rho, args.sigma, args.histogram)

    if args.pdf_type== "G_r":
        df_out = plot_G_r(df, rho, args.sigma, args.histogram)

    if args.output:
        compound = os.path.splitext(args.input_file)[0]
        filename = f"{compound}_{args.pdf_type}.csv"
        df_out.to_csv(filename)
        print(f"Output file saved as {filename}")
        print("")
    
    print("Done\n")


def cli():
    '''Command line interface function'''
    parser = argparse.ArgumentParser(description="xPDFsim")

    # Specification of input file path
    parser.add_argument("input_file", help="cif file")

    # Optional arguments
    parser.add_argument("-o", "--output",
                        action="store_true",
                        help="Optional: if set, pdfs will be written to a csv file in the current directory")

    parser.add_argument("-p", "--pdf_type",
                        type=str,
                        default="g_r",
                        choices=["g_r", "G_r", "R_r"],
                        help="Type of xPDF to be simulated. Options: g_r, G_r, R_r. Default: g_r")
    
    parser.add_argument("-r", "--r_max",
                        type=int,
                        default=20,
                        help="Maximum distance in Angstrom to which the PDF will be calculated. Default: 20")
    
    parser.add_argument("-s", "--sigma",
                        type=float,
                        default=0.1,
                        help="Standard deviation value used for broadening of the PDF histograms. Default: 0.1")
    
    parser.add_argument("-his", "--histogram",
                        action="store_true",
                        help="If set, gausian broadening is disabled and the raw histogram will be plotted/exported.")

    parser.add_argument("-b", "--bin_width",
                        type=float,
                        default=0.05,
                        help="Width of the bins of the PDF histogram in Angstrom. Default: 0.05")


    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
