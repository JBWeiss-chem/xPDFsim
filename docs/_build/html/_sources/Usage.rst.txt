############
Usage
############

You can simply simulate a PDF by opening a terminal in a folder containing your cif file (here *example.cif*) typing:

.. code-block::

    xPDFsim example.cif

Make sure that the conda environment in which you installed xPDFsim is activated. This will generate an .html file of the graph, which is automatically opened with your default browser. On the right, you can use the tools from the tool bar to look at the graph in more detail. You can click on parts of the legend to hide the corresponding plot. The command ``xPDFsim`` has a bunch of options which are listed with:

.. code-block::

  xPDFsim -h
    
One example for using options:

.. code-block::

    xPDFsim -p G_r -o example.cif

This will generate the  :math:`G(r)` function and save it as a csv file.

