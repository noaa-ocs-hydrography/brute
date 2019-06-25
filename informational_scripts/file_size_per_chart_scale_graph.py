# -*- coding: utf-8 -*-

"""
Objective
---------
This generates a graph portraying the tile size relative to chart scale recommened by the IHo verses what NOAA has implemented in the chart rescheme.

Created on Tue April 15 13:44:14 2019

@author: sarah.wolfskehl
"""

import matplotlib.pyplot as plt

# S102 Tile Size Recommendations from Table 11.1 – Informative Grid Resolution and Resulting Tile Size at Chart Scale

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

avgmbpernode = []

# NOAA Band 4 Grid size is 0.3 degrees/18 NM (for scales 1:40,000 and 1:80,000)
# NOAA Band 5 Grid size is 0.075 degrees/4.5 NM (for scales 1:10,000 and 1:20,000)

# y = m * x + b   (m=slope)

m = (10000000 - 256000000) / ((600 * 600) - (5700 * 5700))  # 7.6564 bytes/nodes  m = (y-y)/(x-x)

b = 7243697  # bytes (solved for b by entering y, m, and x for 10MB and 600x600 nodes into line equation)

#  1 nautical mile = 1852 meters

Band4_5m_grid_nodes = ((18 * 1852) / 5) * ((
                                                   18 * 1852) / 5)  # converts nautical miles to meters, divides by grid resolution and mutiplies by itself to get total nodes
Band5_5m_grid_nodes = ((4.5 * 1852) / 5) * ((
                                                    4.5 * 1856) / 5)  # converts nautical miles to meters, divides by grid resolution and mutiplies by itself to get total nodes

Band4_estimated_bytes = m * Band4_5m_grid_nodes + 7243697
Band5_estimated_bytes = m * Band5_5m_grid_nodes + 7243697

plt.plot(10000, Band5_estimated_bytes, 'g^', label="NOAA Band 5 1:10,000")
plt.plot(20000, Band5_estimated_bytes, 'r^', label="NOAA Band 5 1:20,000")
plt.plot(40000, Band4_estimated_bytes, 'b^', label="NOAA Band 4 1:40,000")
plt.plot(80000, Band4_estimated_bytes, 'y^', label="NOAA Band 4 1:80,000")
plt.ylabel('File Size (Bytes)')
plt.xlabel('Scale')
plt.grid()
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.title("Estimated 5m Grid Size for Band 4 and 5")
plt.savefig('graph2.png')
