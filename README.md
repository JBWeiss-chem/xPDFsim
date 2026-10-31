# xPDFsim

Simple simulation of X-ray pair distribution functions. Read the documentation [here](https://xpdfsim.readthedocs.io/en/latest/).

## Features

- Built-in CIF reader (symmetry expansion, partial occupancies, mixed sites, several structures per file)
- Simulation of R(r), ρ(r), g(r), G<sub>K</sub>(r), G(r), D(r), T(r) and N(r), defined as in Peterson *et al.*, [J. Appl. Cryst. 54, 317 (2021)](https://doi.org/10.1107/S1600576720015630)
- X-ray and neutron PDFs
- Simulation of all element-element partial PDFs
- Interactive graphs in the browser
- Export of the simulated xPDFs to a csv file
- Graphical user interface with a 3D structure viewer and instant switching between all of these functions

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install xPDFsim. Its only dependency is numpy.

```bash
pip install xPDFsim
```

## Usage

```bash
xPDFsim example.cif
```

or start the graphical user interface, which opens in your browser:

```bash
xPDFsim-gui
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
