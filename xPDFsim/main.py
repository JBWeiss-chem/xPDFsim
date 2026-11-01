import os
import argparse
import sys
from .cif import read_structures
from .pdf import CONVENTIONS, RADIATIONS, baselines, convert_R_r, gen_R_r, pair_weights, write_csv
from .elements import D_NATURAL
from .utils import broaden_pdfs


def broaden_or_truncate(table, sigma, r_max, histogram, baseline=None):
    """Applies Gaussian broadening, or only cuts the raw histogram at r_max.

    Args:
        table (dictionary): Histograms with a "radii" entry (see pdf.py).
        sigma (float): Standard deviation for Gaussian broadening.
        r_max (float or int): Maximum r to which the PDF was calculated
        histogram (Boolean): If True, no gaussian broadening is applied.
        baseline (function): Optional, takes a grid r and returns the smooth
        part of every function (label -> array), which is not broadened.

    Returns:
        dictionary: Broadened or truncated data.
    """
    if histogram:
        inside = table["radii"] <= r_max
        return {label: values[inside] for label, values in table.items()}
    if baseline is None:
        return broaden_pdfs(table, sigma, r_max)
    peaks = {"radii": table["radii"]}
    peaks.update((label, table[label] - b) for label, b in baseline(table["radii"]).items())
    out = broaden_pdfs(peaks, sigma, r_max)
    for label, b in baseline(out["radii"]).items():
        out[label] = out[label] + b
    return out


def simulate_pdf(structure, convention, r_max, dr, sigma, histogram, radiation="xray", deuteration=D_NATURAL):
    """Calculates a real-space function of a structure.

    Args:
        structure (Structure): Structure of which the PDF shall be calculated.
        convention (string): One of CONVENTIONS.
        r_max (float or int): Maximum r to which the PDF shall be calculated.
        dr (float): Bin width of the histogram.
        sigma (float): Standard deviation for Gaussian broadening.
        histogram (Boolean): If True, no gaussian broadening is applied.
        radiation (string): "xray" or "neutron".
        deuteration (float): Percentage of D among the H atoms (neutrons).

    Returns:
        dictionary: Table of the (broadened) function and its partials.
        float: Number density of the material.
        float: <b_coh>^2 in barn.
    """
    # the Gaussians reach about 6 sigma beyond r_max
    r_extra = 0.0 if histogram else max(1.0, 6 * sigma)
    table, rho = gen_R_r(structure, r_max, dr=dr, radiation=radiation, r_extra=r_extra, deuteration=deuteration)
    weights, mean_b2 = pair_weights(structure, radiation, deuteration)
    converted = convert_R_r(table, rho, convention, weights, mean_b2)
    labels = [label for label in converted if label != "radii"]
    baseline = lambda r: baselines(convention, r, rho, weights, mean_b2, labels)
    return broaden_or_truncate(converted, sigma, r_max, histogram, baseline), rho, mean_b2


def show_in_gui(args, text):
    """Shows the structure and its PDF in the graphical user interface, with the settings of the
    command line, until Ctrl+C is pressed."""
    from .server import serve
    settings = {"pdf_type": args.pdf_type, "radiation": args.radiation, "r_max": args.r_max,
                "dr": args.bin_width, "sigma": args.sigma, "histogram": args.histogram,
                "deuteration": args.deuteration}
    preload = {"file": os.path.basename(args.input_file), "text": text, "settings": settings}
    serve(preload=preload)


def main(args):
    """This function gets executed when the xPDFsim command is run."""
    print(" ")
    try:
        print(f"Input file received: {args.input_file}\n")
        with open(args.input_file, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        structure = read_structures(text)[0][0][1]
    except FileNotFoundError:
        print(f"Error: The file '{args.input_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: The file could not be read as a CIF file: {e}", file=sys.stderr)
        sys.exit(1)

    print("\nSuccessfully processed input file.\n")
    try:
        table, rho, mean_b2 = simulate_pdf(structure, args.pdf_type, args.r_max,
                                           args.bin_width, args.sigma,
                                           args.histogram, args.radiation, args.deuteration)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if args.output:
        compound = os.path.splitext(os.path.basename(args.input_file))[0]
        suffix = "_neutron" if args.radiation == "neutron" else ""
        if args.radiation == "neutron" and args.deuteration != D_NATURAL and "H" in structure.composition():
            suffix += f"_D{args.deuteration:g}"
        filename = f"{compound}_{args.pdf_type}{suffix}.csv"
        write_csv(table, filename)
        print(f"Output file saved as {filename}")
        print(" ")
    if not args.no_plot:
        print("Plotting in the browser...\n")
        show_in_gui(args, text)
    print("Done\n")


def positive_float(value):
    """argparse type that only accepts floats larger than zero."""
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive number")
    return number


def percentage(value):
    """argparse type that only accepts numbers from 0 to 100."""
    number = float(value)
    if not 0 <= number <= 100:
        raise argparse.ArgumentTypeError(f"{value} is not between 0 and 100")
    return number


def cli():
    """Command line interface function."""
    parser = argparse.ArgumentParser(description="xPDFsim: simulation of X-ray and neutron pair distribution functions (PDFs) from CIF files")

    # Specification of input file path
    parser.add_argument("input_file", help="CIF file")

    # Optional arguments
    parser.add_argument("-o", "--output",
                        action="store_true",
                        help="Optional: if set, pdfs will be written to a csv file in the current directory")

    parser.add_argument("-p", "--pdf_type",
                        type=str,
                        default="g_r",
                        choices=list(CONVENTIONS),
                        help="Type of PDF to be simulated, as defined by Peterson et al., J. Appl. Cryst. 54, 317 (2021). Options: G_r, g_r, D_r, GK_r, R_r, rho_r, T_r, N_r. Default: g_r")
    
    parser.add_argument("--radiation",
                        default="xray",
                        choices=list(RADIATIONS),
                        help="Radiation whose scattering lengths weight the pairs: xray (b = Z r_e) or neutron (coherent scattering lengths). Default: xray")

    parser.add_argument("--deuteration",
                        type=percentage,
                        default=D_NATURAL,
                        help=f"Neutrons: percentage of D among the H atoms of the structure (sites given as D stay D). Default: {D_NATURAL:g} (natural abundance)")

    parser.add_argument("-r", "--r_max",
                        type=positive_float,
                        default=20,
                        help="Maximum distance in Angstrom to which the PDF will be calculated. Default: 20")
    
    parser.add_argument("-s", "--sigma",
                        type=positive_float,
                        default=0.1,
                        help="Standard deviation value used for broadening of the PDF histograms. Default: 0.1")
    
    parser.add_argument("-his", "--histogram",
                        action="store_true",
                        help="If set, gaussian broadening is disabled and the raw histogram will be plotted/exported.")

    parser.add_argument("--no-plot",
                        action="store_true",
                        help="Do not show the result in the browser (e.g. for scripts that only write csv files with -o).")

    parser.add_argument("-b", "--bin_width",
                        type=positive_float,
                        default=0.01,
                        help="Rarely worth changing. Width of the bins of the PDF histogram in Angstrom. Default: 0.01")

    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    cli()
