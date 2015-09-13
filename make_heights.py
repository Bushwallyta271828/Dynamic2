# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *

default_stop = 500
default_samples = int(1e6)

def avg_height_find(n, samples=default_samples):
    """
    Author: Xander
    This method computes zeta(n).
    It makes "samples" different arrays of length
    n, where each element in each array is a normally
    distributed variable with mean 0 and standard deviation
    1. It computes the average of max(array) - min(array)
    for each array in our set of "samples" arrays.
    It heavily relies on numpy.
    """
    data = normal(0, 1, (n, samples))
    return average(amax(data, axis=0)
                 - amin(data, axis=0))
    
def make_file(stop=default_stop):
    """
    Author: Xander
    This method calls avg_height_find(...) for every n
    in the range [start, stop], including stop, and writes
    the output to xis.txt in the format:
        n : (xi(n) = n/zeta(n))
    """
    try:
        g = open("xis.txt")
        lines = g.readlines()
        g.close()
        line = lines[-1]
        start = int(line.split(":")[0]) + 1
    except IOError: #xis.txt doesn't exist yet.
        start = 2
    f = open("xis.txt", "a")
    for n in range(start, stop + 1):
        f.write(str(n) + ":" + str(n / avg_height_find(n)) + "\n")
        print n
    f.close()
    
if __name__=="__main__":
    make_file()
