#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 14:57:33 2020

@author: kanthanakorn
"""
import os
import threading

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
# =============================================================================
# Import from TSP
# =============================================================================

home = os.path.expanduser('~')
SExdata = Table.read(os.path.join(home, 'Documents/Data/GrandMa_alert/IM_20200128_212940560_000002_36109501.cat'), format='ascii.sextractor')
CAndata = Table.read(os.path.join(home, 'Documents/Data/GrandMa_alert/IM_20200128_212940560_000002_36109501.candi.dat'), format='ascii')
Cdata = SkyCoord(SExdata['X_WORLD'], SExdata['Y_WORLD'], frame = 'fk5', unit = (u.deg, u.deg))

def Search_Catalog(Cdata_candi):
    confirm_candi1 = Vizier(column_filters={"Rmag":"< 20"},
                            row_limit=-1).query_region(Cdata_candi,
                                                       10*u.arcsec,
                                                       catalog = "I/297",
                                                       verbose=False)
    if str(confirm_candi1) != "Empty TableList":
            USNO_in_radius = SkyCoord(confirm_candi1[0]["RAJ2000"],
                                      confirm_candi1[0]["DEJ2000"],
                                      frame = 'icrs',
                                      unit = (u.deg, u.deg))
            AngSept = USNO_in_radius.separation(Cdata_candi)
    else:
        AngSept = None
    print(AngSept)
    
    return AngSept


if __name__ == '__main__':
    for i in Cdata:
        obj = threading.Thread(target=Search_Catalog, args=(i,))
        obj.start()

