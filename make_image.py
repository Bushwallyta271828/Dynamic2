# -*- coding: utf-8 -*-
from __future__ import division
from pylab import *
from numpy import *
from matplotlib.pyplot import plot, show
from program import *

def create_image():
    lightcurve = normal(0, 1, (999,)) + ((arange(999) - 500)**2) / 100000
    compartmentalization = compartmentalize(lightcurve, max_pval=0.2)[0]
    plot(lightcurve)
    lowers = []
    uppers = []
    xs = []
    print compartmentalization[0]
    for pointnum, point in enumerate(compartmentalization[:-1]):
        next_point = compartmentalization[pointnum + 1]
        block = lightcurve[point: next_point]
        minimum = min(block)
        maximum = max(block)
        xs += [point, next_point - 1]
        lowers += [minimum, minimum]
        uppers += [maximum, maximum]
        if (next_point // 100) > (point // 100):
            plot([100 * (next_point // 100) - 1, 100 * (next_point // 100) - 1,
                  100 * (next_point // 100), 100 * (next_point // 100),
                  100 * (next_point // 100) + 1, 100 * (next_point // 100) + 1],
                 [maximum, minimum, maximum, minimum, maximum, minimum], 'm')
    plot(xs, lowers)
    plot(xs, uppers)
    show()
    
if __name__=="__main__":
    create_image()