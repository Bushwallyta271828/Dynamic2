# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *

def convert():
    """
    Author: Xander
    Given a file heights.txt in the old
    v2.0 format, this function creates the new
    file xis.txt in the new v2.1 format.
    """
    f = open("heights.txt")
    lines = f.readlines()
    f.close()
    g = open("xis.txt", "w")
    for line in lines:
        n = int(line.split(":")[0])
        zeta = float(line.split(":")[1])
        g.write(str(n) + ":" + str(n / zeta) + "\n")
    g.close()

if __name__=="__main__":
    convert()
