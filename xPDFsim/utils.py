import numpy as np


def broaden_pdfs(table, sigma, r_max):
    """Turns histograms into broadened functions.

    Args:
        table (dictionary): Table of histograms with a "radii" entry (see pdf.py).
        sigma (float): Standard deviation for Gaussian broadening.
        r_max (float or int): Maximum r to which the PDF was calculated

    Returns:
        dictionary: Table of the broadened functions, on a grid of 0.02 Å.
    """
    radii = table["radii"]
    dr = radii[1] - radii[0]
    labels = [label for label in table if label != "radii"]
    pdfs = np.column_stack([table[label] for label in labels])
    radii_broad = np.linspace(0, r_max, int(round(r_max * 50)) + 1)

    # Sum of Gaussians centred at every histogram bin, evaluated in chunks to 
    # keep the memory footprint of the kernel matrix small for large r_max
    pdfs_broad = np.empty((len(radii_broad), len(labels)))
    chunk_size = 256
    for start in range(0, len(radii_broad), chunk_size):
        stop = start + chunk_size
        diff = radii_broad[start:stop, np.newaxis] - radii[np.newaxis, :]
        pdfs_broad[start:stop] = np.exp(-diff**2 / (2 * sigma**2)) @ pdfs
    pdfs_broad *= dr / (sigma * np.sqrt(2 * np.pi))

    out = {"radii": radii_broad}
    out.update({label: np.round(pdfs_broad[:, k], 4) for k, label in enumerate(labels)})
    return out
