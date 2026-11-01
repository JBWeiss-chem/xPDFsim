# Changelog

## Unreleased

### New

- Deuteration for the neutron functions: the percentage of D among the H atoms, in the GUI (box 3,
  per structure, applied instantly) and on the command line (`--deuteration`). The default is the
  natural abundance (0.015 %); the field is disabled for structures without hydrogen.

## 0.1.0

### New

- Graphical user interface in the browser (`xPDFsim-gui`): CIF upload and drag & drop, 3D structure
  viewer, interactive PDF plot, export as CSV, SVG and PNG. The command line tool shows its result in it.
- All real-space functions of Peterson *et al.*, J. Appl. Cryst. 54, 317 (2021), with the corrigendum:
  R(r), ρ(r), g(r), G<sub>K</sub>(r), G(r), D(r), T(r) and N(r), with Faber–Ziman weighted partials and
  the asymptotes of Table 4 (`-p`).
- Neutron PDFs with coherent scattering lengths (`--radiation neutron`).
- Highlighting of the atom pairs that contribute most to a region of the PDF in the 3D viewer, with their
  share of the area under R(r); chemically nearly identical distances are grouped.
- Ten example structures from the Crystallography Open Database.
- Partial occupancies and mixed sites are shown in the viewer: the occupied slice of a sphere in the
  atom's colour, the unoccupied rest much lighter.
- `--no-plot` to only calculate, e.g. in scripts with `-o`.
- GitHub workflow that publishes a GitHub release on PyPI (trusted publishing).
- Unit tests (`python -m unittest discover tests`), run by a GitHub workflow on every push.

### Changed

- The command line tool shows the result in the browser GUI instead of a bokeh plot.
- Built-in CIF reader instead of pymatgen; the only dependency is numpy (pymatgen, pandas, bokeh and
  scipy are no longer needed).
- The partials in the csv files are ordered by the first appearance of the elements in the CIF file.
- Only CIF files are read.
- Cell list for the neighbour search and a grid for merging sites in the CIF reader: both scale
  linearly with the size of the structure (e.g. 5000 sites: bonds in 0.02 s instead of 18 s).
- The GUI server only answers to loopback host names and refuses calculations requested by other
  web pages.
- A warning for CIF loops whose number of values is not a multiple of their tags.

### Fixed

- Wrong X-ray normalisation when the order of the elements in the structure differed from the formula.
- CIF files with oxidation states, sites shared by several elements, or files with embedded data (e.g.
  the refinement data that Olex2 writes into `_iucr_refine_fcf_details`) could not be read.
- With more than ten element pairs, the partials beyond the tenth were not plotted.
- `-o` wrote the csv file next to the cif file instead of into the current directory, with an extra
  index column.
- Histogram radii are the bin centres instead of the left bin edges.
- Broadening no longer bends the smooth term of the functions (e.g. −⟨b⟩² of G<sub>K</sub>(r),
  −4πρ₀r of G(r)) near r = 0; only the peaks are broadened.
- The broadened functions are no longer rounded to 4 decimals, which set small partials to zero.
- With wide Gaussians, the functions near r max include the peaks beyond it: the histograms go 6σ
  beyond r max (σ is limited to 0.4 Å in the GUI).
- Cells without a volume (e.g. impossible angles) are rejected instead of giving NaN or hanging.
- Quoted values and text fields that start with `_`, `loop_` or `data_` are read as values.
