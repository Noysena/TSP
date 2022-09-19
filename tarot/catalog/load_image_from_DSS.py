#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 17:06:53 2019

@author: noysena

To download FITS form NASA
"""

import os
import requests

class LOAD_DSS_IMAGE(object):
    def MAIN():
        return("This class will download image from \
               http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl")

    def __init__(self):
        self.home = os.path.expanduser('~')
        self.link = "http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
        self.link_ok = ''
        self.payload = ''


#path = os.path.join(home + "/Documents/web/images/png/7")
#path = os.path.join(home + "/Documents/web")
#link = "http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"

    def Read_Log(self):
        """ Read log on observing night"""
        try:
            self.link_ok = requests.get(self.link, timeout=2)
        except:
            self.link_ok = False
            print("%s --> Error to get link " %self.link)
            pass
#        return self.link_ok

# =============================================================================
# Post paramameter
#Position=%s,%s&Size=%s,%s&Pixels=%s,%s&Rotation=%s&
#Survey=DSS&Scaling=Linear&Projection=Tan&Coordinates=J2000&Return=FITS
# =============================================================================
    def Load_Fits_DSS(self, filename, RA, DE):
        self.payload = {'Position'   : '%s, %s' %(RA, DE),
                'Return'     : 'FITS',
                'Pixels'     : '300,300', #100 less pix, 200 also less, 300 to large size, try 250
                'Size'       : '0.25,0.25',  # FoV in degree 0.2 small, 0.3 lsave_gifarge
                'Rotation'   : '0',
                'Survey'     : 'DSS',
                'Scaling'    : 'Linear',
                'Projection' : 'Tan',
                'Coordinates': 'J2000'
                }
        
        dss_fits = open(filename, 'wb')
        if self.link_ok:
            retrieve_link = requests.post(self.link, params = self.payload)
            dss_fits.write(bytearray(retrieve_link.content))
        dss_fits.close()

# =============================================================================
# Test run
# =============================================================================
if __name__ == '__main__':
    from astropy.table import Table
    home = os.path.expanduser('~')
    p = os.path.join(home, "Documents/Data/GrandMa_alert/IM_20190225_192237528_000001_76279808.candi.dat")
    COOR = Table.read(p, format='ascii')
    name = os.path.basename(p)
    pp = os.path.dirname(p)
    
    test = LOAD_DSS_IMAGE()
    test.Read_Log()

    for i in COOR:
        fitsname = ("%s%s.fits" %(p.replace('.candi.dat', '-dss-'), i["NUMBER"]))
        save_to_file = os.path.join(pp, fitsname)
        test.Load_Fits_DSS(save_to_file, float(i["X_WORLD"]), float(i["Y_WORLD"]))