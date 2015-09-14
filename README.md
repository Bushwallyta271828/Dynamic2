# Dynamic2

This repository contains code for Dynamic Aperture and Dynamic Blocking algorithms as pertaining to Astronomical Photometry; see the more detailed "Dynamic Writeup.pdf" within for an explanation of what these algorithms are.  
There are two folders in this repository: Aperture and Bins. They contain the two algorithms, respectively.
 * Aperture contains the Dynamic Aperture algorithm.
   * aperture.py contains all the code.
   * daperture.py provides a command-line interface for the program.
 * Bins contains the Dynamic Binning algorithm.
   * compartmentalize.py contains all the code for the Binning algorithm.
   * xis.txt is a file containing pre-computed numbers.
   * make_xis.py creates the xis.txt file.
   * dbinning.py provides a command-line interface for the program.
   * convert.py is a program that converts an older file (heights.txt in v2.0) to the newer standard (xis.txt).  
Check the in-code documentation for the fine details.
