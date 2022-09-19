#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 23:01:38 2019

@author: noysena
"""

import os
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.table import Table
from astroquery.imcce import Skybot
from astroquery.xmatch import XMatch
from astropy.coordinates import SkyCoord

class Object_Search(object):
    def __init__(self, input_data):
        self.input_data = input_data
        self.data_cat = Table.read(self.input_data, format='ascii.sextractor')
        self.new = self.input_data.replace('cat', 'new')
        self.fit = self.input_data.replace('cat', 'fits')
        self.xtable = ''
        self.cobj = ''
        self.head = ''
    
    def fits_header(self):
        if os.path.exists(self.new):
            with fits.open(self.new) as hdu:
                self.head = hdu[0].header
        else:
            with fits.open(self.fit) as hdu:
                self.head = hdu[0].hearder
                
    def Xmatch(self, catalog='2MASS'):
        #retrieve catalog
        catalogs = {
        'GaiaDR2' : 'I/345/gaia2',
        'SDSSDR12' : 'V/147/sdss12',
        '2MASS' : 'II/246/out',
        'USNOB1' : 'I/284/out',
        'USNOA2' : 'I/252/out',
        'GLADE2' : 'VII/281/glade2',
        'PanstarrsDR1' : 'II/349/ps1'
        }
        cata = catalogs[catalog]

        #manage data input from sextractor
        data = self.data_cat['X_WORLD', 'Y_WORLD']
        data.write('/tmp/data_input.csv', format='csv', overwrite=True)
        
        self.xtable = XMatch.query(cat1=open('/tmp/data_input.csv'),
                         cat2='vizier:%s' %cata,
                         max_distance=10*u.arcsec,
                         colRA1='X_WORLD',
                         colDec1='Y_WORLD')
        self.xtable.write('/tmp/data_output.csv', format='csv', overwrite=True)
        
    def search_mobj(self, loc, ra, dec):
        """
        field is SkyCoord
        field = SkyCoord(RA*u.deg, DEC*u.deg)
        IAU location {tca : 010, tch : 181, tre 809}
        """
        field = SkyCoord(ra*u.deg, dec*u.deg)
        obs_end_time = Time(self.head['DATE-END']) #Time in 'isot'
        obs_time = Time(obs_end_time.iso, format='iso')
        try:
            self.cobj = Skybot.cone_search(field,
                                           rad = 10*u.arcsec,
                                           epoch = obs_time,
                                           location = loc,
                                           position_error = 120*u.arcsec,
                                           find_planets = True,
                                           find_comets = True,
                                           cache = True)
        except(RuntimeError):
            self.cobj = ("Object Not found: %s %s" %(field.to_string(), obs_end_time))
        
        print(self.cobj)
        

      
if __name__ == '__main__':
    import copy
    home = os.path.expanduser('~')
    input_cat = os.path.join(home, 'Documents/Data/Test/IM_20191207_023356000_000002_26943902.cat')
    obs = Object_Search(input_cat)
    oud = obs.Xmatch(catalog='GaiaDR2')
    ind = obs.data_cat
    obs.fits_header()
    print(obs.xtable)

    #plot angular distance
#        import matplotlib.pyplot as plt
#        plt.figure()
#        plt.hist(obs.xtable['angDist'], bins=100)
#        plt.show()
#       
    test_search = copy.copy(obs.xtable[5670:5674])
    for i in test_search:
        obs.search_mobj('tca', i['ra'], i['dec'])
