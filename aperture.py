# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *
from matplotlib import cm
from scipy.stats import norm
from astropy.io import fits

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
