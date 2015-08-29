# -*- coding: utf-8 -*-
from __future__ import division
import sys
from pylab import *
from numpy import *
from aperture import *

def compute_intensities():
    """
    Author: Xander
    This function is meant to compute star intensities
    given command-line arguments. The arguments that must
    be passed in are:
        1.) the purpose (static aperture? dynamic aperture? both?)
            if purpose is either dynamic or both (i.e. we need a pvalue)
                    |                            |
                    |                            |
                    |                            |
                    |No                          |Yes
                    |                            |
                    |                            |
                    |                            |
            2.) the star's x coordinate        2.) the pvalue
            3.) the star's y coordinate        3.) the star's x coordinate
            4.) the star's bounding radius     4.) the star's y coordinate
            5.) the file to output the results 5.) the star's bounding radius
            6...N) the .fits files             6.) the file to output the results
                                               7...N) the .fits files.
    For either static or dynamic purposes, the output will be of the format:
        intensity[0](newline)
        intensity[1](newline)
        intensity[3](newline)
        ...
    For "both", the output will be of the format:
        dynamic_intensity[0] static_intensity[0](newline)
        dynamic_intensity[1] static_intensity[1](newline)
        dynamic_intensity[2] static_intensity[2](newline)
        ...
    where the dynamic and static intensities
    are separated by single spaces. 
    """
    args = sys.argv
    purpose = args[1]
    if purpose in ["dynamic", "both"]:
        pval = float(args[2])
        x = int(args[3])
        y = int(args[4])
        rad = int(args[5])
        output = args[6]
        places = args[7: len(args)]
    else: #purpose = "static"
        x = int(args[2])
        y = int(args[3])
        rad = int(args[4])
        output = args[5]
        places = args[6: len(args)]
    if purpose == "static":
        intensities = []
        for place in places:
            photo = Photo(place)
            photo.load()
            intensities.append(photo.old_intensity((x, y), rad))
            photo.close()
        f = open(output, "w")
        for intensity in intensities:
            f.write(str(intensity) + "\n")
        f.close()
    elif purpose == "dynamic":
        intensities = []
        for place in places:
            photo = Photo(place)
            photo.load()
            intensities.append(photo.new_intensity((x, y), rad, pval))
            photo.close()
        f = open(output, "w")
        for intensity in intensities:
            f.write(str(intensity) + "\n")
        f.close()
    else: #purpose = "both"
        new_intensities = []
        old_intensities = []
        for place in places:
            photo = Photo(place)
            photo.load()
            new_intensities.append(photo.new_intensity((x, y), rad, pval))
            old_intensities.append(photo.old_intensity((x, y), rad))
            photo.close()
        f = open(output, "w")
        for intensitynum in range(len(new_intensities)):
            f.write(str(new_intensities[intensitynum])
                  + " "
                  + str(old_intensities[intensitynum])
                  + "\n")
        f.close()
        
if __name__=="__main__":
    compute_intensities()