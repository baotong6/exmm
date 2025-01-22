'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-09-30 15:50:55
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-09-30 17:14:16
FilePath: /code/gaia_match.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import sys
import os
from astropy import units as u
from astropy.coordinates import SkyCoord,match_coordinates_sky
from astroquery.gaia import Gaia
from astropy.coordinates import Angle
from astropy.table import Table
from astroquery.vizier import Vizier
from matches_info import output_xmmmatchesinfo
# Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
# print(Gaia)
(matches,xmmmatchinfolist)=output_xmmmatchesinfo()
def matches_togaia():
    ra_e=matches['RA'];dec_e=matches['DEC']
    coord_e = SkyCoord(ra_e*u.degree, dec_e*u.degree)
    maxrad = 0.3* u.degree                      
    result = Vizier(row_limit = -1,columns=['*','Date']).query_region(coord_e[0], radius=maxrad, catalog="I/355")
    gaia_table = result[0]
    coord_gaia = SkyCoord(gaia_table['RAJ2000'],gaia_table['DEJ2000'])
    # Define max sep in arcsec to consider good match
    distance_limit = 5 * u.arcsec
    for i in range(len(coord_e)):
        idxc, idxcatalog, d2d, d3d = coord_gaia.search_around_sky(coord_e[i], distance_limit)
        print(coord_e[idxc][0],coord_gaia[idxcatalog],d2d.arcsec)
        
if __name__ == '__main__':
    matches_togaia()


