# Dynamic2
Code for Dynamic Aperture and Dynamic Blocking algorithms as pertaining to Astronomical Photometry; see the more detailed "Dynamic Writeup.pdf" within for an explanation of what these algorithms are.  
The \__doc\__ strings should be completely self-explanatory in figuring out how the code behaves on a programming level.  
program.py is the file with the code in it. 
All pictures must be .fits files.
heights.txt is a list of pre-computed values stored as a .txt. They store what I called the "zeta" function in the .pdf file. DO NOT STORE heights.txt IN A SEPERATE DIRECTORY TO program.py AND EXPECT compartmentalize() TO WORK!!! (Unless you MANUALLY specifiy barmap as a keyword argument)  
daperture.py and dbinning.py are interfaces that don't actually compute anything - they just pass the task on to program.py. They are useful for interfacing with the code, even if you don't know Python or want to figure out which functions should be called. Their \__doc\__ strings should be self-explanatory regarding how to format command-line arguments and how to expect the output file to be formatted.
