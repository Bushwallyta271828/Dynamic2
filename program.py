from __future__ import division
from pylab import *
from matplotlib import cm
from scipy.stats import norm
from astropy.io import fits

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

def box((x, y), rad, shape, border=0):
    """
    Author: Xander
    box is used by the Photo class constantly.
    Given a point, a desired radius, and the shape
    of the array of pixels, box provides the box 
    around (x, y) that doesn't go closer than border
    pixels to the edge. This is used to protect stars
    close to the edge of the frame from causing Index out
    of Bounds exceptions and the like.
    """
    x_start = max(x - rad, border)
    x_stop = min(x + rad, shape[0] - border)
    y_start = max(y - rad, border)
    y_stop = min(y + rad, shape[1] - border)
    return ((x_start, x_stop), (y_start, y_stop))
    
def recurse(point, forbidden, points_mask):
    """
    Author: Xander
    recurse(...) is used by the Photo class for Dynamic Aperture, 
    but since it doesn't use or call anything about the Photo class,
    nor does it have any real connection on its own to astronomy,
    I didn't put it in with the other methods.
    Given a point on an abstract grid, a mask of places not to go,
    and a mask of points allowed to visit by the function,
    it returns a tuple containing:
        1.) a list of the points that can be reached by going up,
            down, left, and right without leaving the points_mask
            region or entering the forbidden region.
        2.) A new forbidden region, consisting of the old forbidden
            region and those points now visited.
    """
    neighbors = [(point[0] - 1, point[1]),
                 (point[0] + 1, point[1]),
                 (point[0], point[1] - 1),
                 (point[0], point[1] + 1)]
    forbidden.add(point)
    region = [point]
    for neighbor in neighbors:
        if (neighbor not in forbidden
            and 0 <= neighbor[0] <= len(points_mask) - 1
            and 0 <= neighbor[1] <= len(points_mask[0]) - 1
            and points_mask[neighbor]):
            recurse_result = recurse(neighbor, forbidden, points_mask)
            region += recurse_result[0]
            forbidden = recurse_result[1]
    return (region, forbidden)

class Photo:
    def __init__(self, name):
        """
        Author: Xander
        This class represents a .fits photograph taken through
        a telescope of the sky. It is passed the name of the file
        where this photograph exists.
        """
        self.name = name

    def load(self):
        """
        Author: Xander
        This method loads the photograph into RAM.
        As well as creating self.headers and self.hdulist,
        it makes the self.scidata numpy array of the pixels in the 
        image, which is used constantly by the program.
        """
        self.hdulist = fits.open(self.name)
        self.scidata = self.hdulist[0].data
        self.scidata = self.scidata.transpose()
        if len(self.scidata.shape) == 3:
            blank_axis = list(self.scidata.shape).index(1)
            self.scidata = self.scidata.sum(axis=blank_axis)
        self.headers = self.hdulist[0].header

    def close(self):
        """
        Author: Xander
        close(self) deletes the arrays created 
        in load() to free up RAM. This is needed when
        many .fits files are opened and closed in the duration
        of one run of the program.
        """
        try: self.hdulist.close()
        except: pass
        try: del self.scidata
        except: pass
        try: del self.hdulist
        except: pass
        try: del self.headers
        except: pass

    def stretch_plot(self, lower, upper):
        """
        Author: Xander
        Given an upper limit to be set to 255 
        and a lower limit to be set to 0, this method plots
        a stretched version of self.scidata.
        """
        new_scidata = self.scidata.copy()
        new_scidata[new_scidata > upper] = upper
        new_scidata[new_scidata < lower] = lower
        new_scidata = (new_scidata - lower) / (upper - lower)
        imshow(new_scidata.transpose(), cmap = cm.Greys_r)

    def auto_plot(self):
        """
        Author: Xander
        This method finds default values of a lower limit
        and an upper limit based on self.scidata, then hands
        the plotting task over to stretch_plot(...).
        """
        data_pool = self.scidata.reshape(len(self.scidata)*len(self.scidata[0]))
        data_pool.sort()
        lower = data_pool[int(len(data_pool) * 0.00025)]
        upper = data_pool[int(len(data_pool) * 0.99975)]
        self.stretch_plot(lower, upper)
        
    def new_intensity(self, (x, y), rad, pval=0.05):
        """
        Author: Xander
        new_intensity is the method that will be interfaced with
        by the user for Dynamic Aperture. Given a point (x, y) and a
        radius rad in which a star resides (Note: the star isn't of
        radius rad, but rather the star is completely contained in 
        a circle of radius rad and center (x, y)), this method
        finds the intensity of the star. pval is an optional argument that
        is relevant to the find_star_mask(...) call. The appropriate (Dynamic)
        fraction of the average background noise is taken out before returning
        the intensity.
        """
        ((x_start, x_stop), (y_start, y_stop)) = box((x, y), rad, self.scidata.shape)
        xpos, ypos = meshgrid(arange(x_start, x_stop), arange(y_start, y_stop))
        xpos = xpos.transpose()
        ypos = ypos.transpose()
        star_mask = self.find_star_mask((x, y), rad, xpos, ypos, pval)
        not_star_mask = logical_not(star_mask)
        inner_size = sum(star_mask)
        outer_size = sum(not_star_mask)
        inner_sum = sum(self.scidata[xpos, ypos] * star_mask)
        outer_sum = sum(self.scidata[xpos, ypos] * not_star_mask)
        return max(inner_sum - outer_sum * inner_size / outer_size, 0)
        
    def old_intensity(self, (x, y), rad):
        """
        Author: Xander
        old_intensity computes the intensity in a non-dynamic way.
        It assumes the star is somewhere within a circle of center (x, y)
        and radius rad, and computes the intensity with the backgroud taken
        as the region between a square around the circle and the circle itself.
        """
        ((x_start, x_stop), (y_start, y_stop)) = box((x, y), rad, self.scidata.shape)
        xpos, ypos = meshgrid(arange(x_start, x_stop), arange(y_start, y_stop))
        xpos = xpos.transpose()
        ypos = ypos.transpose()
        star_mask = ((xpos - x)**2 + (ypos - y)**2 < rad**2)
        not_star_mask = logical_not(star_mask)
        inner_size = sum(star_mask)
        outer_size = sum(not_star_mask)
        inner_sum = sum(self.scidata[xpos, ypos] * star_mask)
        outer_sum = sum(self.scidata[xpos, ypos] * not_star_mask)
        return max(inner_sum - outer_sum * inner_size / outer_size, 0)

    def outer_stats(self, (x, y), rad):
        """
        Author: Xander
        outer_stats is used to get a representative feel
        of the background levels in a particular patch of the image.
        Given a circle of center (x, y) and radius rad, outer_stats 
        finds the mean and the standard deviation of pixels bewteen
        a square of center (x, y) and sidelength 2*rad and the specified
        circle.
        """
        ((x_start, x_stop), (y_start, y_stop)) = box((x, y), rad, self.scidata.shape)
        xvals, yvals = meshgrid(arange(x_start - x, x_stop - x),
                                arange(y_start - y, y_stop - y))
        xvals = xvals.transpose()
        yvals = yvals.transpose()
        mask = (xvals**2 + yvals**2 > rad**2)
        scidata_slice = (self.scidata[x_start: x_stop, y_start: y_stop])[mask]
        size = sum(mask)
        mean = sum(scidata_slice) / size
        std = sqrt(sum((scidata_slice - mean)**2) / (size - 1))
        return (mean, std)
        
    def cutoff(self, (x, y), rad, pval=0.05):
        """
        Author: Xander
        Given a p-value and a circle centered at (x, y)
        with radius rad, cutoff(...) finds the cutoff
        such that pval fraction of the background data points would
        be above that cutoff, if we assume that background points 
        have pixel intensities distributed as a gaussian.
        """
        mean, std = self.outer_stats((x, y), rad)
        zscore = norm.ppf(1 - pval)
        return mean + zscore*std
    
    def find_star_mask(self, (x, y), rad, xpos, ypos, pval=0.02):
        """
        Author: Xander
        find_star_mask is the backbone of the Dynamic Aperture program.
        Given a circle with center (x, y), radius rad, x-values xpos,
        y-values ypos, and a P-value, find_star_mask finds what 
        we call the "star" in the Dynamic Aperture program. It does 
        this by finding all points above a cutoff specified by the
        cutoff(...) program, and then finding the largest one
        connected by moving up, down, right, and left on the grid.
        """
        cutoff = self.cutoff((x, y), rad, pval)
        original_mask = (self.scidata[xpos, ypos] >= cutoff)
        xposes = (xpos[original_mask] - xpos[0, 0])
        yposes = (ypos[original_mask] - ypos[0, 0])
        points = zip(xposes, yposes)
        best = []
        forbidden = set()
        for point in points:
            if point not in forbidden:
                region, forbidden = recurse(point,
                                            forbidden,
                                            original_mask)
                if len(region) > len(best):
                    best = region
        if best == []: return zeros(xpos.shape, dtype=int)
        else:
            allowed_mask = zeros(xpos.shape, dtype=int)
            for p in best:
                allowed_mask[p] = 1
            return allowed_mask

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