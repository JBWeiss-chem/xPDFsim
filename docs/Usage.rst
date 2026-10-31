############
Usage
############

You can simply simulate a PDF by opening a terminal in a folder containing your cif file (here *example.cif*) typing:

.. code-block::

    xPDFsim example.cif

This calculates the PDF and shows it, together with the structure, in the graphical user interface (see below), which opens in your default browser; press Ctrl+C in the terminal to stop it. Use ``--no-plot`` to only calculate (e.g. in scripts together with ``-o``). The command ``xPDFsim`` has a bunch of options which are listed with:

.. code-block::

  xPDFsim -h
    
One example for using options:

.. code-block::

    xPDFsim -p G_r -o example.cif

This will generate the  :math:`G(r)` function and save it as a csv file. Add ``--radiation neutron`` for the neutron PDF.


Graphical user interface
========================

xPDFsim also comes with a graphical user interface that runs in your browser. Start it with:

.. code-block::

    xPDFsim-gui

Open one or more cif files with the *Open CIF…* button or drop them onto the window, or pick one of ten classic inorganic, organic and hybrid structures from the Crystallography Open Database in the *Examples…* menu. The structure is shown in a 3D viewer, and its PDF is calculated right away. The panel on the right switches between all functions of the Introduction, changes the Gaussian broadening and selects the partial PDFs without recalculating anything. Only a new maximum distance or bin width requires a recalculation. Drag across the plot to zoom in, double-click to show the full range, and export the curves as csv, svg or png files. Left-click the plot to highlight the atom–atom distances of the region around that distance in the 3D viewer, each with its share of the area under :math:`R(r)` in that region. The highlighting stays while you turn the structure; click into the highlighted region or press Esc to remove it. The width of the region and the minimum share of a highlighted distance are set in the PDF options. The supercell shown is set in the toolbar of the viewer; *≈ r max* makes it about as thick as the maximum distance of the PDF.
