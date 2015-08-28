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
                     x=0,
                     memo=None,
                     max_length=float("inf"),
                     barmap=None,
                     max_pval=0.2,
                     min_size=None):
    """
    Author: Xander
    compartmentalize(...) is the backbone of the Dynamic Bins program.
    Given a lightcurve, a starting frame number, a memo table, the maximum
    number of frames that can be put in a bin, list of (index=number
    of items in a bin) to (value_at_index=expectation value of maximum -
    minimum given a standard deviation of 1), the maximum fraction of the time
    a measurement may lie outside the error bars, and the minimum number of 
    frames a bin can contain, this function returns a tuple containing:
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
            and not the standard deviations themselves, we can use barmap to
            infer the standard deviation (assuming the lightcurve is distributed
            as a gaussian for small intervals) with the formula:
                standard_deviation = height / barmap[number of frames]
            Since we must weigh the bins proportionally with the number of frames
            they contain, we arrive at the badness formula:
                badness of a bin = (number of frames) * height / barmap[number of frames]
            The returned partitioning is the partitioning, given that 
            no bin can have more than max_length items in it, that minimizes the total 
            badness.
        3.) The new memo table. This third returned item is only helpful to further
            recursive calls of the compartmentalize function - it is of no direct use to
            the user.
    Default values mean that the user only HAS to input a lightcurve. If not set,
    barmap will be pulled heights.txt, which has this exact list written out as a 
    .txt file for reasonable bin sizes (starting with 2). The only reasonable optional
    keywords for a user to specify are max_length and max_pval.
    The program uses Dynamic Programming to find the best partitioning in O(n) time,
    where n is the length of the lightcurve. Since any algorithm must at least look
    at every intensity of the lightcurve, this program provides the best
    concieveable big-O. However, the constant is quite large.
    Another note: The algorithm is only O(n) because of the max_length 
    requirement, which saves us from looking at blocks stretching very long
    distances. If max_length is set to float("inf"), for instance, the algorithm
    becomes O(n^2). Think of the program as O(n * min(n, max_length - min_size)).
    The last two keyword arguments passed in are linked. The user specifies max_pval,
    the maximum probability that a measurement can land outside the error bars, which
    the computer then uses to find min_size, the minimum number of frames allowed in
    a bin. Since there is no need to re-compute this quantity, it is passed in directly
    to recursive calls. The formula for min_size in terms of max_pval is:
        min_size = max(2, int(2 / max_pval))
    This can be derived by some combination of math and quick computer scripts to 
    check that math for small values of min_size.
    """
    if memo == None: memo = {}
    if barmap == None:
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
    if min_size == None:
        min_size = max(2, int(2 / max_pval))
        """
        Inverting the formula pval = 2 / (min_size + 1), but since int
        rounds down, the -1 must be discarded.
        """
    if x in memo: return memo[x] + (memo,)
    elif len(lightcurve) - x < min_size:
        partition = [x, len(lightcurve)]
        badness = float("inf")
        return (partition, badness, memo)
    elif len(lightcurve) - x in [2, 3]:
        partition = [x, len(lightcurve)]
        block = lightcurve[x:]
        height = max(block) - min(block)
        badness = len(block) * height / barmap[len(block)]
        return (partition, badness, memo)
    else:
        pre_lightcurve_slice = lightcurve[x: x + min_size]
        maximum = max(pre_lightcurve_slice)
        minimum = min(pre_lightcurve_slice)
        best_badness = float("inf")
        best_partition = []
        for new_x in range(x + min_size, min(len(lightcurve) - 1, x + max_length + 1)):
            height = maximum - minimum
            new_badness = (new_x - x) * height / barmap[new_x - x]
            compartmentalization = compartmentalize(lightcurve,
                                                    new_x,
                                                    memo,
                                                    max_length,
                                                    barmap,
                                                    max_pval,
                                                    min_size)
            memo = compartmentalization[2]
            found_badness = compartmentalization[1]
            badness = new_badness + found_badness
            if badness < best_badness:
                best_badness = badness
                best_partition = [x] + compartmentalization[0]
            if lightcurve[new_x] < minimum:
                minimum = lightcurve[new_x]
            if lightcurve[new_x] > maximum:
                maximum = lightcurve[new_x]
        if len(lightcurve) - x <= max_length:
            maximum = max(maximum, lightcurve[-1])
            minimum = min(minimum, lightcurve[-1])
            height = maximum - minimum
            new_badness = (len(lightcurve) - x) * height / barmap[len(lightcurve) - x]
            if new_badness < best_badness:
                best_badness = new_badness
                best_partition = [x, len(lightcurve)]
        memo[x] = (best_partition, best_badness)
        return (best_partition, best_badness, memo)