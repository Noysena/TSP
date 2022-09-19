#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:55:03 2019

@author: noysena
"""

import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
from astropy.table import QTable
from astropy.table import vstack
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord

class OBJECT_SEARCH(object):
    def __init__(self):
        self.tbdata = ''
        self.cobj = ''
        self.head = ''
        self.loc = ''
        self.new = ''
        self.fits = ''
        self.obs_time=''
        
    def Asteroid_Comet(coor, radius, loc, obs_time):
        try:
            cobj = Skybot.cone_search(coor,
                                           rad = radius*u.deg,
                                           epoch = obs_time,
                                           location = loc,
                                           position_error = 120*u.arcsec,
                                           find_planets = True,
                                           find_comets = True,
                                           find_asteroids = True,
                                           cache = False)
        except:
            #Handling exceptions
            cobj = sys.exc_info()
            
        return cobj
        
    def Check_Asteroid_Comet(self, file_candi_dat):
        """
        input data is *.candi.dat (astropy table ASCII format)
        field is SkyCoord
        field = SkyCoord(RA*u.deg, DEC*u.deg)
        
        change timeout from 300 to 3000
        /home/noysena/.local/lib/python3.6/site-packages/astroquery/imcce
        """
        
        self.tbdata = Table.read(file_candi_dat, format='ascii')
        self.new = file_candi_dat.replace('.candi.dat', '.new')
        self.fits = file_candi_dat.replace('.candi.dat', '.fits')
        
        if os.path.exists(self.new):
            with fits.open(self.new) as hdu:
                self.head = hdu[0].header
        else:
            with fits.open(self.fit) as hdu:
                self.head = hdu[0].header
                
        if self.head['TELESCOP'] == 'TAROT CALERN':
            self.loc = '010'
        elif self.head['TELESCOP'] == 'TAROT CHILI':
            self.loc = '181'
        elif self.head['TELESCOP'] == 'TAROT REUNION':
            self.loc = '809'
        else:
            print("Not TAROT telescope")
            sys.exit()
                
        field = SkyCoord(self.tbdata['X_WORLD']*u.deg, self.tbdata['Y_WORLD']*u.deg)
        
        for coor in field:
            obs_end_time = Time(self.head['DATE-END']) #Time in 'isot'
            self.obs_time = Time(obs_end_time.iso, format='iso')
            print(self.obs_time)
        
            try:
                self.cobj = Skybot.cone_search(coor,
                                               rad = 0.088*u.deg,
                                               epoch = self.obs_time,
                                               location = self.loc,
                                               position_error = 120*u.arcsec,
                                               find_planets = True,
                                               find_comets = True,
                                               find_asteroids = True,
                                               cache = False)
                print(self.cobj)
            except:
                #Handling exceptions
                self.cobj = sys.exc_info()

    def Check_Asteroid_Comet_FoV(self, fitsfile):
            with fits.open(fitsfile) as hdu:
                head = hdu[0].header
                wcs = WCS(hdu[0].header)

            RA = head['RA']
            DEC = head['DEC']
            
            field = SkyCoord(RA*u.deg, DEC*u.deg)
            obs_end_time = Time(head['DATE-END']) #Time in 'isot'
            obs_time = Time(obs_end_time.iso, format='iso')
            print(obs_time)
            try:
                head['TELESCOP']
            except(KeyError):
                head['TELESCOP'] = 'TAROT REUNION'
                
            if head['TELESCOP'] == 'TAROT CALERN':
                radius = 1.29
                loc = '010'
                print("%s %s" %(head['TELESCOP'], field.to_string()))
                cobj = OBJECT_SEARCH.Asteroid_Comet(field, radius, loc, obs_time)
                
            elif head['TELESCOP'] == 'TAROT CHILI':
                radius = 1.29
                loc = '181'
                print("%s %s" %(head['TELESCOP'], field.to_string()))
                cobj = OBJECT_SEARCH.Asteroid_Comet(field, radius, loc, obs_time)

            elif head['TELESCOP'] == 'TAROT REUNION':
                radius = 1.29
                loc = '809'

            #Divided to 4 part which equals to 2k but vstack in the end
                c1 = wcs.all_pix2world(1024, 1024, 0)
                field1 = SkyCoord(c1[0]*u.deg, c1[1]*u.deg)
                cobj1 = OBJECT_SEARCH.Asteroid_Comet(field1, radius, loc, obs_time)
                if type(cobj1) is QTable:
                    pass
                else:
                    cobj1 = ''
                
                c2 = wcs.all_pix2world(1024, 3072, 0)
                field2 = SkyCoord(c2[0]*u.deg, c2[1]*u.deg)
                cobj2 = OBJECT_SEARCH.Asteroid_Comet(field2, radius, loc, obs_time)
                if type(cobj2) is QTable:
                    pass
                else:
                    cobj2 = ''

                c3 = wcs.all_pix2world(3072, 3072, 0)
                field3 = SkyCoord(c3[0]*u.deg, c3[1]*u.deg)
                cobj3 = OBJECT_SEARCH.Asteroid_Comet(field3, radius, loc, obs_time)
                if type(cobj3) is QTable:
                    pass
                else:
                    cobj3 = ''
                
                c4 = wcs.all_pix2world(3072, 1024, 0)
                field4 = SkyCoord(c4[0]*u.deg, c4[1]*u.deg)
                cobj4 = OBJECT_SEARCH.Asteroid_Comet(field4, radius, loc, obs_time)
                if type(cobj4) is QTable:
                    pass
                else:
                    cobj4 = ''
                cobj = [cobj1, cobj2, cobj3, cobj4]
                
                for i in cobj:
                    if isinstance(i, str):
                        cobj.remove(i)
                    else:
                        pass
                    
                cobj = vstack(cobj)
            else:
                print("Not TAROT telescope")
                sys.exit()
                
            return cobj
                

if __name__ == '__main__':
    home = os.path.expanduser('~')
    filename = os.path.join(home, "Documents/Data/Test/IM_20191207_023356000_000002_26943902.candi.dat")
    fitsfile = filename.replace('.candi.dat', '.fits')
    run = OBJECT_SEARCH()
    cobj = run.Check_Asteroid_Comet_FoV(fitsfile)
    print(cobj)
