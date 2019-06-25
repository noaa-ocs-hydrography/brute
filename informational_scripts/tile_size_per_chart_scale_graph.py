# -*- coding: utf-8 -*-

"""
Objective
---------
This generates a graph portraying the tile size relative to chart scale recommened by the IHo verses what NOAA has implemented in the chart rescheme.

Created on Tue April 15 13:44:14 2019

@author: sarah.wolfskehl
"""

import matplotlib.pyplot as plt

grid_resolution_list = [900,
                        900,
                        450,
                        210,
                        105,
                        54,
                        27,
                        13,
                        6,
                        3,
                        2,
                        1,
                        1,
                        1,
                        1]
scale = [10000000,
         3500000,
         1500000,
         700000,
         350000,
         180000,
         90000,
         45000,
         22000,
         12000,
         8000,
         4000,
         3000,
         2000,
         1000]

resulting_tile_size_at_10MB = [291,
                               291,
                               145,
                               68,
                               34,
                               17.5,
                               8.7,
                               4.2,
                               1.9,
                               1.0,
                               0.6,
                               0.3,
                               0.3,
                               0.3,
                               0.3]

resulting_tile_size_at_256MB = [2770,
                                2770,
                                1385,
                                646,
                                323,
                                166,
                                83,
                                40,
                                18.5,
                                9.0,
                                6.0,
                                3.0,
                                3.0,
                                3.0,
                                3.0]

plt.plot(10000, 4.5, 'r^', label="NOAA Band 5 1:10,000 ")
plt.plot(20000, 4.5, 'g^', label="NOAA Band 5 1:20,000")
plt.plot(40000, 18, 'b^', label="NOAA Band 4 1:40,000")
plt.plot(80000, 18, 'c^', label="NOAA Band 4 1:80,000")
plt.plot(160000, 72, 'm^', label="NOAA Band 3 1:160,000")
plt.plot(320000, 72, 'y^', label="NOAA Band 3 1:320,000")
plt.plot(scale, resulting_tile_size_at_10MB, 'g', label="IHO Recommended Minimum")
plt.plot(scale, resulting_tile_size_at_256MB, 'b', label="IHO Recommended Maximum")
plt.ylabel('Tile Size (nautical miles)')
plt.xlabel('Scale')
plt.xlim((5000, 330000))
plt.ylim((0, 80))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title("RESCHEME BAND SIZE PER SCALE")
plt.savefig('graph.png')
