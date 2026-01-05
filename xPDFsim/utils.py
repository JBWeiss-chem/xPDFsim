import numpy as np
import pandas as pd
import re
from bokeh.models import TeX


def get_fractions(formula):
    """
    Counts the total number of elements and calculates the fraction of each element
    in a chemical formula given as a string.
    """
    formula = formula.replace(' ', '')

    # Regex to find element symbols and their counts
    pattern = r'([A-Z][a-∏z]*)(\d*\.?\d*)'
    elements = re.findall(pattern, formula)

    element_counts = {}
    total_atoms = 0.0
    for element, count in elements:
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



def broaden_pdfs(df, sigma, r_max):
    """
    Turns histograms into broadened functions. The broadening increases with r.
    """
    radii = df["radii"].to_numpy()
    dr = radii[1] - radii[0]
    labels = df.columns.to_list()
    labels.remove('radii')
    radii_broad = np.linspace(0, r_max, r_max*50)

    # Put everything into a dataframe
    data = {'radii': radii_broad}
    df_out = pd.DataFrame(data)

    for label in labels:
        pdf = df[label].to_numpy()
        pdf_broad = np.zeros_like(radii_broad)
        for r, y in zip(radii, pdf):
            norm_factor = dr / (sigma * np.sqrt(2 * np.pi))
            pdf_broad += y * norm_factor * np.exp(-((radii_broad - r)**2) / (2 * sigma**2))
        df_out[label] = pdf_broad
    df_out = df_out[df_out['radii'] <= r_max]
    return df_out.round(4)


def format_plot(p):
    """
    Formats the Brokeh plot nicely.
    """
    p.xaxis.axis_label_text_font_size = "20pt"
    p.yaxis.axis_label_text_font_size = "20pt"
    p.yaxis.major_label_text_font_size = "20pt"
    p.xaxis.major_label_text_font_size = "20pt"
    p.legend.label_text_font_size = "20pt"
    p.legend.glyph_width = 70
    p.legend.click_policy="hide"
    p.xaxis.axis_label = TeX(r"r / \mathrm{\mathring{A}}")
