# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *
from numpy import *
from compartmentalize import *
from aperture import *
from matplotlib import cm
from scipy.stats import norm
from astropy.io import fits

def get_lightcurves((x, y) = (289, 42), rad=20):
    """
    Author: Xander
    This function calls Photo.new_intensity and
    Photo.old_intensity for every frame in a dataset
    located at:
        C:\\Users\\Nicholas\\Desktop\\1718+48\\<framenumber>.fit
    where <framenumber> is preceeded by the appropriate number of 
    0s to make the total length of the digits between the \\ and the
    .fit be 7. The default star under examination is located at 
    (x=289, y=42) and is completely encompassed by a circle at that 
    location of radius rad=20. Those parameters can be tweaked with
    keywords to focus in on another star if needed.
    """
    new_values = []
    old_values = []
    for i in range(1, 1713):
        photo = Photo("C:\\Users\\Nicholas\\Desktop\\1718+48\\"
                    + "0"*(7 - len(str(i))) + str(i) + ".fit")
        photo.load()
        new_values.append(photo.new_intensity((x, y), rad, pval=0.02))
        old_values.append(photo.old_intensity((x, y), rad))
        photo.close()
    return (new_values, old_values)
        

def plot_gen(dist,
             include_list=[True, False, False, True],
             style_list=["solid", "solid", "solid", "solid"]):
    """
    Author: Xander
    plot_gen is the specific program I use to make my graph.
    It does not implement any "algorithms" or do anything
    to reduce error bars - it only calls methods and plots results.
    It doesn't return anything.
    It takes as arguments a "dist" parameter - the maximum
    number of frames that can be grouped together with Dynamic Blocking,
    and an include_list parameter, which specifies which of:
        1.) Dynamic Aperture and Dynamic Bins
        2.) Dynamic Aperture but No Dynamic Bins
        3.) No Dynamic Aperture but Dynamic Bins
        4.) No Dynamic Aperture and No Dynamic Bins
    to graph. style_list represents what style each of these graphs
    is to be done in.
    """
    new_values, old_values = get_lightcurves()
    if include_list[0]:
        compartmentalization = compartmentalize(new_values)
        xs = []
        ys = []
        for xnum, x in enumerate(compartmentalization[0][:-1]):
            y = compartmentalization[0][xnum + 1]
            block = new_values[x: y]
            minimum = min(block)
            maximum = max(block)
            height = maximum - minimum
            xs += [x, y - 1]
            ys += [height, height]
        plot(xs,
             ys,
             color="b",
             linestyle=style_list[0],
             label="Dynamic Aperture and Dynamic Bins")
    if include_list[1]:
        xs = []
        ys = []
        for x in range(0, 1712, dist):
            if x + dist <= 1712:
                y = x + dist
                block = new_values[x: y]
                minimum = min(block)
                maximum = max(block)
                height = maximum - minimum
                xs += [x, y - 1]
                ys += [height, height]
        plot(xs,
             ys,
             color="g",
             linestyle=style_list[1],
             label="Dynamic Aperture but No Dynamic Bins")
    if include_list[2]:
        compartmentalization = compartmentalize(old_values)
        xs = []
        ys = []
        for xnum, x in enumerate(compartmentalization[0][:-1]):
            y = compartmentalization[0][xnum + 1]
            block = old_values[x: y]
            minimum = min(block)
            maximum = max(block)
            height = maximum - minimum
            xs += [x, y - 1]
            ys += [height, height]
        plot(xs,
             ys,
             color="r",
             linestyle=style_list[2],
             label="No Dynamic Aperture but Dynamic Bins")
    if include_list[3]:
        xs = []
        ys = []
        for x in range(0, 1712, dist):
            if x + dist <= 1712:
                y = x + dist
                block = old_values[x: y]
                minimum = min(block)
                maximum = max(block)
                height = maximum - minimum
                xs += [x, y - 1]
                ys += [height, height]
        plot(xs,
             ys,
             color="m",
             linestyle=style_list[3],
             label="No Dynamic Aperture and No Dynamic Bins")
    ylim(ymin=0)
    legend(fontsize=18)
    xlabel("Frame Number", fontsize=24, labelpad=30)
    ylabel("Error Bar Height", fontsize=24, labelpad=40)
    suptitle("Error Bars for Dynamic Aperture and Dynamic Bins (Static Bin Size = "
           + str(dist)
           + ")",
             fontsize=36)
    show()
    
def lightcurve_plot():
    """
    Author: Xander
    This function creates a graph of the lightcurve, where 
    the data before 1712 // 2 is shown with both Dynamic Blocking
    and Dynamic Aperture, while the data after is shown with neither.
    It is meant to demonstrate the stark contrast that both these methods 
    bring to the error bars of the data set, side by side.
    """
    new_values, old_values = get_lightcurves()
    plot(range(1712), new_values[:1712 // 2] + old_values[1712 // 2:], "b", label="Lightcurve")
    plot([1712 // 2,
          1712 // 2],
         [max(max(new_values),
              max(old_values)),
          min(min(new_values),
              min(old_values))], "m", label="Division Line")
    partition = compartmentalize(new_values[:1712 // 2])[0]
    xs = []
    uppers = []
    lowers = []
    for placenum, place in enumerate(partition[:-1]):
        next_place = partition[placenum + 1]
        maximum = max(new_values[place: next_place])
        minimum = min(new_values[place: next_place])
        xs += [place, next_place - 1]
        uppers += [maximum, maximum]
        lowers += [minimum, minimum]
    for i in range(1712 // 2, 1712, 20):
        if i + 20 < 1712:
            maximum = max(old_values[i: i + 20])
            minimum = min(old_values[i: i + 20])
            xs += [i, i + 19]
            uppers += [maximum, maximum]
            lowers += [minimum, minimum]
    plot(xs, uppers, "r", label="Upper Bound")
    plot(xs, lowers, "g", label="Lower Bound")
    legend(fontsize=18)
    xlabel("Frame Number", fontsize=24, labelpad=30)
    ylabel("Intensity", fontsize=24, labelpad=40)
    suptitle("Lightcurve with and without both Dynamic Aperture and Dynamic Bins",
             fontsize=36)
    show()
    
def average_length():
    """
    Author: Xander
    average_length finds and returns the average length of 
    the Dynamic Bins algorithm, when set to find the BEST
    partitioning, which happens when max_length = float("inf")
    (i.e. compartmentalize is O(n^2)).
    """
    new_values = get_lightcurves()[0]
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    average_length = 0
    size = 0
    for partition_point_num, partition_point in enumerate(partitioning[:-1]):
        next_partition_point = partitioning[partition_point_num + 1]
        length = next_partition_point - partition_point
        average_length += length
        size += 1
    average_length /= size
    return average_length
    
def average_weighted_length():
    """
    Author: Xander
    average_length finds and returns the average length of 
    the Dynamic Bins algorithm, when set to find the BEST
    partitioning, which happens when max_length = float("inf")
    (i.e. compartmentalize is O(n^2)).
    """
    new_values = get_lightcurves()[0]
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    average_length = 0
    for partition_point_num, partition_point in enumerate(partitioning[:-1]):
        next_partition_point = partitioning[partition_point_num + 1]
        length = next_partition_point - partition_point
        average_length += length * length
    average_length /= len(new_values)
    return average_length
    
def five_number_summary():
    """
    Author: Xander
    This function finds the five number summary (minimum, first quartile,
    median, third quartile, and maximum) of the lenghts of bins returned by
    the Dynamic Binning algorihthm above.
    """
    new_values = get_lightcurves()[0]
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    lengths = []
    for partition_point_num, partition_point in enumerate(partitioning[:-1]):
        next_partition_point = partitioning[partition_point_num + 1]
        length = next_partition_point - partition_point
        lengths.append(length)
    lengths.sort()
    l = len(lengths)
    minimum = lengths[0]
    first_quartile = lengths[l // 4]
    median = lengths[l // 2]
    third_quartile = lengths[3 * l // 4]
    maximum = lengths[-1]
    return (minimum, first_quartile, median, third_quartile, maximum)

def five_number_weighted_summary():
    """
    Author: Xander
    This function finds the five number summary (minimum, first quartile,
    median, third quartile, and maximum) of the lenghts of bins, each
    repeated the number of frames they hold, returned by the Dynamic
    Binning algorihthm above.
    """
    new_values = get_lightcurves()[0]
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    lengths = []
    for partition_point_num, partition_point in enumerate(partitioning[:-1]):
        next_partition_point = partitioning[partition_point_num + 1]
        length = next_partition_point - partition_point
        lengths += [length]*length
    lengths.sort()
    l = len(lengths)
    minimum = lengths[0]
    first_quartile = lengths[l // 4]
    median = lengths[l // 2]
    third_quartile = lengths[3 * l // 4]
    maximum = lengths[-1]
    return (minimum, first_quartile, median, third_quartile, maximum)
    
def weighted_average_drop_multiple(dist=38):
    """
    Author: Xander
    This function finds, weighted by number of frames in each bin,
    the average decrease in error bar height due to employing 
    Dynamic Binning.
    """
    new_values = get_lightcurves()[0]
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    nothing_sum = 0
    everything_sum = 0
    for x in range(0, 1712, dist):
        if x + dist <= 1712:
            y = x + dist
            block = new_values[x: y]
            height = max(block) - min(block)
            nothing_sum += height * dist
    for pointnum, point in enumerate(partitioning[:-1]):
        next_point = partitioning[pointnum + 1]
        block = new_values[point: next_point]
        height = max(block) - min(block)
        everything_sum += height * (next_point - point)
    return nothing_sum / everything_sum
    
def weighted_average_total_drop_multiple(dist=38):
    """
    Author: Xander
    This function finds, weighted by number of frames in each bin,
    the average decrease in error bar height due to employing 
    Dynamic Binning AND Dynamic Aperture.
    """
    new_values, old_values = get_lightcurves()
    partitioning = compartmentalize(new_values, max_length = float("inf"))[0]
    nothing_sum = 0
    everything_sum = 0
    for x in range(0, 1712, dist):
        if x + dist <= 1712:
            y = x + dist
            block = old_values[x: y]
            height = max(block) - min(block)
            nothing_sum += height * dist
    for pointnum, point in enumerate(partitioning[:-1]):
        next_point = partitioning[pointnum + 1]
        block = new_values[point: next_point]
        height = max(block) - min(block)
        everything_sum += height * (next_point - point)
    return nothing_sum / everything_sum
    
def weighted_average_aperture_drop_multiple(dist=38):
    """
    Author: Xander
    This function finds, weighted by number of frames in each bin,
    the average decrease in error bar height due to employing 
    Dynamic Aperture.
    """
    new_values, old_values = get_lightcurves()
    nothing_sum = 0
    everything_sum = 0
    for x in range(0, 1712, dist):
        if x + dist <= 1712:
            y = x + dist
            block = old_values[x: y]
            height = max(block) - min(block)
            nothing_sum += height * dist
    for x in range(0, 1712, dist):
        if x + dist <= 1712:
            y = x + dist
            block = new_values[x: y]
            height = max(block) - min(block)
            everything_sum += height * dist
    return nothing_sum / everything_sum
    
def badness_drop_factor(dist=38):
    """
    Author: Xander
    This function finds the increase factor by using Static
    Binning to the "badness".
    """
    new_values = get_lightcurves()[0]
    dynamic_badness = compartmentalize(new_values)[1]
    static_badness = 0
    zeta = 4.27932121178
    for x in range(0, 1712, dist):
        if x + dist <= 1712:
            y = x + dist
            block = new_values[x: y]
            height = max(block) - min(block)
            static_badness += dist * height / zeta
    return static_badness / dynamic_badness