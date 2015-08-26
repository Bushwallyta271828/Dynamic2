# -*- coding: utf-8 -*-
from __future__ import division
import sys
from pylab import *
from program import *

def compute_binning():
    """
    Author: Xander
    """
    args = sys.argv
    input_file = args[1]
    output_file = args[2]
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