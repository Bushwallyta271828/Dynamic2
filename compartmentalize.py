# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *
from numpy import *

"""
Comments:
A Block is a Bin.
heights.txt must be kept in the same directory at all times!
"""

def compartmentalize(lightcurve,
                     max_length=float("inf"),
                     max_pval=0.2):
    """
    Author: Xander
    compartmentalize(...) is the backbone of the Dynamic Bins program.
    Given a lightcurve, the maximum number of frames that can be put in a bin,
    and the maximum fraction of the time a measurement may lie
    outside the error bars, this function returns a tuple containing:
        1.) A list of the optimal breaking points.
            In other words, if the first item in the returned tuple reads
            [a, b, c, ..., z], then the best bins to use are [a, b), [b, c),
            ..., [y, z), where "z" must be the length of the lightcurve, and "a"
            must be the starting position given (the keyword "x")
        2.) The total badness of this combination. The badness of a collection
            of bins is defined to be the sum of the badnesses of each bin.
            I define the goal of this program to be maximizing the information
            we know about the data set from the error bars, which in turn I 
            define to mean minimizing the average standard deviation, with larger
            bins proportionally contributing more weight than smaller bins.
            Since a "badness" metric can be multiplied by a constant to give an 
            equally valid "badness" metric giving the exact same results, and
            since len(lightcurve) - x = the number of data points we are breaking
            up does not depend on the way we break the frames up, we can ignore the
            divided by part of the average and just make our problem the minimization
            of the sum of the standard deviations, weighted by size, of the bins. 
            However, since the end result only depends on the error bar heights,
            and not the standard deviations themselves, we can use zeta in heights.txt
            to infer the standard deviation (assuming the lightcurve is distributed
            as a gaussian for small intervals) with the formula:
                standard_deviation = height / barmap[number of frames]
            Since we must weigh the bins proportionally with the number of frames
            they contain, we arrive at the badness formula:
                badness of a bin = (number of frames) * height / barmap[number of frames]
            The returned partitioning is the partitioning, given that 
            no bin can have more than max_length items in it, that minimizes the total 
            badness. For more information on the algorithm, check out Dynamic Writeup.pdf
    The function uses Dynamic Programming to find the best partitioning in O(n) time,
    where n is the length of the lightcurve. Since any algorithm must at least look
    at every intensity of the lightcurve, this program provides the best
    concieveable big-O. However, the constant is quite large.
    Another note: The algorithm is only O(n) because of the max_length 
    requirement, which saves us from looking at blocks stretching very long
    distances. If max_length is set to float("inf"), for instance, the algorithm
    becomes O(n^2). Think of the program as O(n * min(n, max_length - min_size)).
    """
    max_length = min(max_length, len(lightcurve))
    f = open("heights.txt")
    lines = f.readlines()
    f.close()
    barmap = [None, None] #Neither 0 nor 1 may be called - they have no meaning.
    for line in lines:
        height = float(line[:-1].split(":")[1])
        barmap.append(height)
    slope = (barmap[-1] - barmap[-21]) / 20
    while max_length > len(barmap) - 1:
        barmap.append(barmap[-1] + slope)
    min_size = max(2, int(2 / max_pval))
    memovalues = [None]*len(lightcurve) + [len(lightcurve)]
    memovalues[len(lightcurve) - min_size + 1:-1] = [([], float("inf"))]*min_size
    i = len(lightcurve) - min_size
    while i >= 0:
        minimum = min(lightcurve[i:i + min_size])
        maximum = max(lightcurve[i:i + min_size])
        best_badness = float("inf")
        best_partitioning = []
        for j in range(i + min_size, min(len(lightcurve), i + max_length + 1)):
            new_badness = (j - i) * (maximum - minimum) / barmap[j - i]
            found_partitioning, found_badness = memovalues[j]
            if new_badness + found_badness > best_badness:
                best_badness = new_badness + found_badness
                best_partitioning = [i] + found_partitioning
            if lightcurve[j] > maximum:
                maximum = lightcurve[j]
            elif lightcurve[j] < minimum:
                minimum = lightcurve[j]
        memovalues[i] = (best_partitioning, best_badness)
    return memovalues[0]