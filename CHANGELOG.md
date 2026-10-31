# Changelog

## 0.1.0

### New

- Graphical user interface in the browser (`xPDFsim-gui`): CIF upload and drag & drop, 3D structure
  viewer, interactive PDF plot, export as CSV, SVG and PNG. The command line tool shows its result in it.
- All real-space functions of Peterson *et al.*, J. Appl. Cryst. 54, 317 (2021), with the corrigendum:
  R(r), ρ(r), g(r), G<sub>K</sub>(r), G(r), D(r), T(r) and N(r), with Faber–Ziman weighted partials and
  the asymptotes of Table 4 (`-p`).
- Neutron PDFs with coherent scattering lengths (`--radiation neutron`).
- Highlighting of the atom pairs that contribute most to a region of the PDF in the 3D viewer, with their
  share of the area under R(r).
- Ten example structures from the Crystallography Open Database.
- Partial occupancies and mixed sites are shown in the viewer: the occupied slice of a sphere in the
  atom's colour, the unoccupied rest much lighter.
- `--no-plot` to only calculate, e.g. in scripts with `-o`.

### Changed

- The command line tool shows the result in the browser GUI instead of a bokeh plot.
- Built-in CIF reader instead of pymatgen; the only dependency is numpy (pymatgen, pandas, bokeh and
  scipy are no longer needed).
- The partials in the csv files are ordered by the first appearance of the elements in the CIF file.
- Only CIF files are read.

### Fixed

- Wrong X-ray normalisation when the order of the elements in the structure differed from the formula.
- CIF files with oxidation states, sites shared by several elements, or files with embedded data (e.g.
  the refinement data that Olex2 writes into `_iucr_refine_fcf_details`) could not be read.
- With more than ten element pairs, the partials beyond the tenth were not plotted.
- `-o` wrote the csv file next to the cif file instead of into the current directory, with an extra
  index column.
- Histogram radii are the bin centres instead of the left bin edges.
