from __future__ import division
from pylab import *

default_start = 5
default_stop = 100

def avg_height_find(n, samples=1000000):
    data = normal(0, 1, (n, samples))
    maxs = amax(data, axis=0)
    mins = amin(data, axis=0)
    heights = maxs - mins
    return average(heights)
    
def make_file(start=default_start, stop=default_stop):
    avg_heights = {}
    f = open("heights.txt", "w")
    for i in avg_heights:
        f.write(str(i) + ":" + str(avg_heights[i]) + "\n")
    for n in range(start, stop + 1):
        if n not in avg_heights:
            f.write(str(n) + ":" + str(avg_height_find(n)) + "\n")
            print n
    f.close()
    
if __name__=="__main__":
    make_file(start=2, stop=50)