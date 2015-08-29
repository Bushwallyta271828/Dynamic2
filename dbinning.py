# -*- coding: utf-8 -*-
from __future__ import division
import sys
from pylab import *
from numpy import *
from compartmentalize import *

def compute_binning():
    """
    Author: Xander
    This function is a command-line interface
    to the compartmentalize(...) function in
    program.py. The command-line arguments
    that must be passed are:
        1.) the input file name
        2.) the output file name
        3.) the maximum acceptible length (inf for infinity)
        4.) the maximum acceptible P-value.
    The output file will be of the format:
        badness(newline)
        compartmentalize_output[0][0](newline)
        compartmentalize_output[0][1](newline)
        compartmentalize_output[0][2](newline)
        ...
    """
    args = sys.argv
    input_file = args[1]
    output_file = args[2]
    if args[3] == "inf":
        max_length = float("inf")
    else:
        max_length = int(args[3])
    max_pval = float(args[4])
    f = open(input_file)
    lines = f.readlines()
    f.close()
    values = []
    for line in lines:
        values.append(float(line[:-1]))
    compartmentalization = compartmentalize(values,
                                            max_length=max_length,
                                            max_pval=max_pval)
    partitioning = compartmentalization[0]
    badness = compartmentalization[1]
    g = open(output_file, "w")
    g.write(str(badness) + "\n")
    for partition in partitioning:
        g.write(str(partition) + "\n")
    g.close()
    
if __name__=="__main__":
    compute_binning()