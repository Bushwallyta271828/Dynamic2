# Dynamic2
This repository contains code for Dynamic Aperture and Dynamic Blocking algorithms as pertaining to Astronomical Photometry; see the more detailed "Dynamic Writeup.pdf" within for an explanation of what these algorithms are.  
For a file-by-file explanation of what each script does from a managment perspective:  
    The Photographs folder contains photographs used in "Dynamic Writeup.pdf". They are not used as any part of the program.  
    .gitattributes and .gitignore are files relating to git and github, not to the algorithms.  
    aperture.py contains all the functions required to run the Dynamic Aperture program.  
    compartmentalize.py contains all the functions required to run the Dynamic Binning program.  
    daperture.py is a command-line interface with the aperture.py program that doesn't require a knowledge of python.  
    dbinning.py is a command-line interface with the compartmentalize.py program that doesn't require a knowledge of python.  
    Dyanmic Writeup.pdf is a file describing the algorithms implemented here.  
    heights.txt is a .txt file of pre-computed quantities used in compartmentalize.py (zeta(n) in the pdf for n in [2, 500])  
    make_heights.py computes heights.txt  
    README.md is this very file! As you can probably guess, it discusses what this repository is and how to use it.  
For a how-to-use-the-code level tutorial, read the \__doc\__ strings attached to each and every function. 
