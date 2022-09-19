#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 16:36:49 2018

@author: kanthanakorn
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord


def Iterate(Cdata, Ccata):
    """
    Input : data coordinate, catalog coordinate; both in ICRS standard
    Output: indices, angular distance, distance(in case coordinate includeed)
            matching of data, data's coordinate
    """
    A, B = 1, 0
    c = 0
    idx, d2d, d3d = Cdata.match_to_catalog_sky(Ccata)
    matches = Ccata[idx]
#    print("Matches indexs\n", matches)
    
    while (A - B) > 0 and c < 10:
        A = np.round(np.median(d2d.arcsec), 5)
        delta_RA = matches.ra.deg - Cdata.ra.deg
        delta_DEC = matches.dec.deg - Cdata.dec.deg
        reRA = Cdata.ra.deg + np.median(delta_RA)
        reDEC = Cdata.dec.deg + np.median(delta_DEC)
        Cdata = SkyCoord(reRA, reDEC, frame = 'icrs', unit = (u.deg, u.deg))
        idx, d2d, d3d = Cdata.match_to_catalog_sky(Ccata)
        matches = Ccata[idx]
        B = np.round(np.median(d2d.arcsec), 5)
        c += 1
        
        try:
            print("Median Sept: %s" %B)
        except(KeyboardInterrupt, SystemExit):
            raise
               
    return idx, d2d, d3d, matches, Cdata

if __name__ == "__main__":
    pass