# Dynamic2
Code for Dynamic Aperture and Dynamic Blocking algorithms as pertaining to Astronomical Photometry; see the more detailed "Dynamic Writeup.pdf" within for an explanation of what these algorithms are.
The \__doc\__ strings should be completely self-explanatory in figuring out how the code behaves on a programming level.
v8.py is the file with the code in it.
heights.txt is a list of pre-computed values stored as a .txt. They store what I called the "zeta" function in the .pdf file. DO NOT STORE heights.txt IN A SEPERATE DIRECTORY TO v8.py AND EXPECT compartmentalize() TO WORK!!! (Unless you MANUALLY specifiy barmap as a keyword argument)
