# xPDFsim

Simulation of X-ray and neutron pair distribution functions (PDFs) from CIF files, with an interactive
graphical user interface in the browser. Read the documentation [here](https://xpdfsim.readthedocs.io/en/latest/).

![The 3D viewer of xPDFsim with the Cu–O bonds of HKUST-1 highlighted](https://raw.githubusercontent.com/JBWeiss-chem/xPDFsim/main/docs/gui_viewer.png)

## Features

- Built-in CIF reader (symmetry expansion, partial occupancies, mixed sites, several structures per file)
- Simulation of R(r), ρ(r), g(r), G<sub>K</sub>(r), G(r), D(r), T(r) and N(r), defined as in Peterson *et al.*, [J. Appl. Cryst. 54, 317 (2021)](https://doi.org/10.1107/S1600576720015630)
- X-ray and neutron PDFs, with all element-element partial PDFs
- Graphical user interface: 3D structure viewer, instant switching between all functions, broadening
  or raw histogram, zoom and pan, export as CSV, SVG and PNG
- Click a region of the PDF to highlight the atom pairs that contribute most to it in the 3D viewer,
  with their share of the PDF
- Ten example structures from the Crystallography Open Database
- Export of the simulated PDFs to a CSV file from the command line
- Only one dependency: numpy

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install xPDFsim. Its only dependency is numpy.

```bash
pip install xPDFsim
```

## Usage

Simulate the PDF of a structure; the result opens in the graphical user interface in your browser:

```bash
xPDFsim example.cif
```

Calculate the neutron G(r) and save it as a CSV file, without opening the browser:

```bash
xPDFsim -p G_r --radiation neutron -o --no-plot example.cif
```

Start the graphical user interface without a file:

```bash
xPDFsim-gui
```

All options are listed with `xPDFsim -h`.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/). The bundled three.js is licensed under the MIT
licence; the example structures from the COD are in the public domain.
