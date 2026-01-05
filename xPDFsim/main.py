import os
import numpy as np
from bokeh.plotting import figure, show
from bokeh.palettes import Category10_10
from bokeh.models import TeX
import argparse
import sys
import pandas as pd
from pymatgen.core import Structure
from .pdf import *
from .utils import *


def plot_R_r(df, rho, sigma, r_max, histogram):
    """
    Turns dataframe into R(r) plot.
    """
    # Apply gaussian broadening
    df_broad = broaden_pdfs(df, sigma, r_max)
    if histogram == True:
        df_broad = df
    radii = df_broad["radii"].to_numpy()
    labels = df.columns.to_numpy()[:-2]

    R_r_total = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, sizing_mode="stretch_both")
    p.yaxis.axis_label = TeX(r"R(r) / \mathrm{\mathring{A}}^{-1}")
    
    # Plot full PDF
    p.line(radii, R_r_total, color="black", line_width=1.5, legend_label="R(r)")

    # # Plot trend line to check validity
    # y = 4*np.pi*rho*radii**2
    # p.line(radii, y, color="red", line_width=1.5, legend_label="check")

    for label, color in zip(labels, Category10_10):
        R_r = df_broad[label].to_numpy()
        p.line(radii, R_r, color=color,line_width=1.5, legend_label=label)

    format_plot(p)
    p.legend.location = "top_left"

    show(p)
    return df_broad


def plot_g_r(df, rho, sigma, r_max, histogram):
    """
    Turns dataframe into g(r) plot.
    """
    radii = df["radii"].to_numpy()
    labels = df.columns.to_numpy()[:-2]

    # Calculate g_r partials from R_r partials
    for label in labels:
        df[label] = df[label] / (4*np.pi*rho*radii**2)
    df['total'] = df['total'] / (4*np.pi*rho*radii**2)

    # Apply gaussian broadening
    df_broad = broaden_pdfs(df, sigma, r_max)
    if histogram == True:
        df_broad = df
    radii_broad = df_broad["radii"].to_numpy()
    g_r_total = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, sizing_mode="stretch_both")
    p.yaxis.axis_label = TeX(r"g(r) / \mathrm{\mathring{A}}^{-2}")

    # Plot full PDF
    p.line(radii_broad, g_r_total, color="black", line_width=1.5, legend_label="g(r)")    

    for label, color in zip(labels, Category10_10):
        g_r = df_broad[label].to_numpy()
        p.line(radii_broad, g_r, color=color,line_width=1.5, legend_label=label)

    p.legend.location = "top_right"
    format_plot(p)
        
    show(p)
    return df_broad
    

def plot_G_r(df, rho, sigma, r_max, histogram):
    """
    Turns dataframe into G(r) plot.
    """
    radii = df["radii"].to_numpy()
    total_R_r = df["total"].to_numpy()

    G_r = total_R_r / radii - 4 * np.pi * rho * radii
    data = {'radii': radii}
    data["total"] = G_r
    df_G_r = pd.DataFrame(data)
    
    # Apply gaussian broadening
    df_broad = broaden_pdfs(df_G_r, sigma, r_max)
    if histogram == True:
        df_broad = df_G_r
    radii_broad = df_broad["radii"].to_numpy()
    G_r_broad = df_broad["total"].to_numpy()

    # Create plot window
    p = figure(title=None, sizing_mode="stretch_both")
    p.yaxis.axis_label = TeX(r"G(r) / \mathrm{\mathring{A}}^{-2}")

    # Plot full PDF
    p.line(radii_broad, G_r_broad, color="black", line_width=1.5, legend_label="G(r)")

    p.legend.location = "top_right"
    format_plot(p)
    show(p)
    return df_broad


def main(args):
    print(" ")
    try:
        print(f"Input file received: {args.input_file}\n")
        structure = Structure.from_file(args.input_file, primitive=True)
    except FileNotFoundError:
        print(f"Error: The file '{args.input_file}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
    
    print("\nSuccessfully processed input file.\n")
    df, rho = gen_R_r(structure, args.r_max, dr=args.bin_width)

    print("Plotting...\n")
    if args.pdf_type == "R_r":
        df_out = plot_R_r(df, rho, args.sigma, args.r_max, args.histogram)

    if args.pdf_type == "g_r":
        df_out = plot_g_r(df, rho, args.sigma, args.r_max, args.histogram)

    if args.pdf_type== "G_r":
        df_out = plot_G_r(df, rho, args.sigma, args.r_max, args.histogram)

    if args.output:
        compound = os.path.splitext(args.input_file)[0]
        filename = f"{compound}_{args.pdf_type}.csv"
        df_out.to_csv(filename)
        print(f"Output file saved as {filename}")
        print(" ")
    print("Done\n")


def cli():
    '''
    Command line interface function.
    '''
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
                        default=0.01,
                        help="Rarely worth changing. Width of the bins of the PDF histogram in Angstrom. Default: 0.01")

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
