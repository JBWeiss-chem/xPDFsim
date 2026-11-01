"""Element data: atomic numbers, coherent neutron scattering lengths, covalent radii
(Cordero et al., Dalton Trans. 2008) and Jmol colours (radii and colours taken from NetsUnhinged)."""

# element symbols in the order of the atomic numbers (D: deuterium, Z = 1)
SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
    "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
    "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
    "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
ATOMIC_NUMBERS = {sym: z for z, sym in enumerate(SYMBOLS, start=1)}
ATOMIC_NUMBERS["D"] = 1

# coherent neutron scattering lengths in fm of the natural elements (Neutron Data Booklet,
# 2nd ed., 2003; Sears' values; as tabulated in pymatgen, MIT licence) and of deuterium
NEUTRON_B = {
    "H": -3.739, "He": 3.26, "Li": -1.9, "Be": 7.79, "B": 5.3, "C": 6.646, "N": 9.36, "O": 5.803,
    "F": 5.654, "Ne": 4.566, "Na": 3.63, "Mg": 5.375, "Al": 3.449, "Si": 4.1491, "P": 5.13, "S": 2.847,
    "Cl": 9.577, "Ar": 1.909, "K": 3.67, "Ca": 4.7, "Sc": 12.29, "Ti": -3.438, "V": -0.3824, "Cr": 3.635,
    "Mn": -3.73, "Fe": 9.45, "Co": 2.49, "Ni": 10.3, "Cu": 7.718, "Zn": 5.68, "Ga": 7.288, "Ge": 8.185,
    "As": 6.58, "Se": 7.97, "Br": 6.795, "Kr": 7.81, "Rb": 7.09, "Sr": 7.02, "Y": 7.75, "Zr": 7.16,
    "Nb": 7.054, "Mo": 6.715, "Tc": 6.8, "Ru": 7.03, "Rh": 5.88, "Pd": 5.91, "Ag": 5.922, "Cd": 4.87,
    "In": 4.065, "Sn": 6.225, "Sb": 5.57, "Te": 5.8, "I": 5.28, "Xe": 4.92, "Cs": 5.42, "Ba": 5.07,
    "La": 8.24, "Ce": 4.84, "Pr": 4.58, "Nd": 7.69, "Pm": 12.6, "Sm": 0.8, "Eu": 7.22, "Gd": 6.5,
    "Tb": 7.38, "Dy": 16.9, "Ho": 8.01, "Er": 7.79, "Tm": 7.07, "Yb": 12.43, "Lu": 7.21, "Hf": 7.7,
    "Ta": 6.91, "W": 4.86, "Re": 9.2, "Os": 10.7, "Ir": 10.6, "Pt": 9.6, "Au": 7.63, "Hg": 12.692,
    "Tl": 8.776, "Pb": 9.405, "Bi": 8.532, "Ra": 10, "Th": 10.31, "Pa": 9.1, "U": 8.417, "Np": 10.55,
    "Am": 8.3,
    "D": 6.671,
}

# natural abundance of D in hydrogen in % (Neutron Data Booklet), the default deuteration
D_NATURAL = 0.015


def hydrogen_b(deuteration=D_NATURAL):
    """Coherent neutron scattering length in fm of hydrogen with the given percentage of D among
    its atoms, interpolated linearly between natural H and D (0 %: 1H, -3.7406 fm)."""
    if not 0 <= deuteration <= 100:
        raise ValueError("The deuteration must be between 0 and 100 %.")
    b_h, b_d = NEUTRON_B["H"], NEUTRON_B["D"]
    return b_h + (deuteration - D_NATURAL) / (100 - D_NATURAL) * (b_d - b_h)


COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42,
    "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15,
    "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98,
    "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87,
    "Lu": 1.87, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21, "Ac": 2.15, "Th": 2.06, "Pa": 2.00,
    "U": 1.96, "Np": 1.90, "Pu": 1.87, "Am": 1.80, "Cm": 1.69, "D": 0.31,
}

NONMETALS = {"H", "D", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar",
             "As", "Se", "Br", "Kr", "Te", "I", "Xe", "At", "Rn"}


def is_metal(el: str) -> bool:
    return el not in NONMETALS


JMOL_COLORS = {
    "H": "#FFFFFF", "D": "#FFFFC0", "He": "#D9FFFF", "Li": "#CC80FF", "Be": "#C2FF00",
    "B": "#FFB5B5", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D", "F": "#90E050",
    "Na": "#AB5CF2", "Mg": "#8AFF00", "Al": "#BFA6A6", "Si": "#F0C8A0", "P": "#FF8000",
    "S": "#FFFF30", "Cl": "#1FF01F", "K": "#8F40D4", "Ca": "#3DFF00", "Sc": "#E6E6E6",
    "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7", "Mn": "#9C7AC7", "Fe": "#E06633",
    "Co": "#F090A0", "Ni": "#50D050", "Cu": "#C88033", "Zn": "#7D80B0", "Ga": "#C28F8F",
    "Ge": "#668F8F", "As": "#BD80E3", "Se": "#FFA100", "Br": "#A62929", "Rb": "#702EB0",
    "Sr": "#00FF00", "Y": "#94FFFF", "Zr": "#94E0E0", "Nb": "#73C2C9", "Mo": "#54B5B5",
    "Ru": "#248F8F", "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F",
    "In": "#A67573", "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00", "I": "#940094",
    "Cs": "#57178F", "Ba": "#00C900", "La": "#70D4FF", "Ce": "#FFFFC7", "Nd": "#C7FFC7",
    "Eu": "#61FFC7", "Gd": "#45FFC7", "Tb": "#30FFC7", "Dy": "#1FFFC7", "Er": "#00E675",
    "Yb": "#00BF38", "Hf": "#4DC2FF", "W": "#2194D6", "Pt": "#D0D0E0", "Au": "#FFD123",
    "Hg": "#B8B8D0", "Pb": "#575961", "Bi": "#9E4FB5", "U": "#008FFF",
}


def color(el: str) -> str:
    return JMOL_COLORS.get(el, "#FF1493")
