#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 15:43:01 2019

@author: noysena
"""

import os
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy import units as u
from astropy.table import vstack
from astropy.table import QTable
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
from tarot.utility.limit_in_fits import Limit_Coor_to_Fits

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
    
def Asteroid_Comet_2(RA, DEC, telescope, obs_time):
    if telescope == 'TAROT CALERN':
        loc = '010'
    elif telescope == 'TAROT CHILI':
        loc = '181'
    elif telescope == 'TAROT REUNION':
        loc = '809'
    else:
        print("Not TAROT telescope")
        sys.exit()
        
    field = SkyCoord(RA*u.deg, DEC*u.deg)
        
    try:
        cobj = Skybot.cone_search(field,
                                       rad = 16*u.arcsec,
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
#    Build condition to get stat that No solar object in circle is True
    try:
        cobj_txt = cobj[1].args[0][0:15]
    except:
        cobj_txt = ''
    if cobj_txt == 'No solar system':
        sobj = True
    else:
        sobj = False
    
    if sobj:
        pass
    else:
        try:
            cobj_txt = cobj[1].args[0]
            print(cobj_txt)
           # print(cobj[1].args[0])
        except(IndexError):
            pass
        
    return sobj
    
def Check_Asteroid_Comet(fitsfile, RA, DEC):
    """
    input data is *.candi.dat (astropy table ASCII format)
    field is SkyCoord
    field = SkyCoord(RA*u.deg, DEC*u.deg)
    
    change timeout from 300 to 3000
    /home/noysena/.local/lib/python3.6/site-packages/astroquery/imcce
    """
    if isinstance(RA, float) and isinstance(DEC, float):
        pass
    elif isinstance(RA, int) and isinstance(DEC, int):
        pass
    elif isinstance(RA, str) and isinstance(DEC, str):
        RA = float(RA)
        DEC = float(RA)
    else:
        print("Not correct RA and DEC %s %s" %(RA, DEC))
        
    with fits.open(fitsfile) as hdu:
        head = hdu[0].header
        
    field = SkyCoord(RA*u.deg, DEC*u.deg)
    
    obs_end_time = Time(head['DATE-END']) #Time in 'isot'
    obs_time = Time(obs_end_time.iso, format='iso')
    if head['TELESCOP'] == 'TAROT CALERN':
        loc = '010'
    elif head['TELESCOP'] == 'TAROT CHILI':
        loc = '181'
    elif head['TELESCOP'] == 'TAROT REUNION':
        loc = '809'
    else:
        print("Not TAROT telescope")
        sys.exit()
        
    
    print("Observatory: %s\nObserved time: %s" %(loc, obs_time))
    print("Coordinate: %s" %field.to_string())
    try:
        cobj = Skybot.cone_search(field,
                                  rad = 0.088*u.deg,
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
#        print(cobj[1])
    
    if type(cobj) is QTable:
        cobj = Limit_Coor_to_Fits(fitsfile, catalog=cobj, ratext='RA', detext='DEC')
    else:
        pass
    
    return cobj

#check full fits with radius 1.29, 4k image use "check_solar_object.py"
def Check_Asteroid_Comet_Full_2k_FITS(fitsfile, RA, DEC, recode=True):
    """
    RA and DEC is center FoV
    input data is *.candi.dat (astropy table ASCII format)
    field is SkyCoord
    field = SkyCoord(RA*u.deg, DEC*u.deg)
    """
    if isinstance(RA, float) and isinstance(DEC, float):
        pass
    elif isinstance(RA, int) and isinstance(DEC, int):
        pass
    elif isinstance(RA, str) and isinstance(DEC, str):
        RA = float(RA)
        DEC = float(RA)
    else:
        print("Not correct RA and DEC %s %s" %(RA, DEC))
    
    with fits.open(fitsfile) as hdu:
        head = hdu[0].header
    
    #       save solar object to votable
    if fitsfile.split('.')[-1] == 'new':
        save_votable = fitsfile.replace('.new', '-sol.vot')
    elif fitsfile.split('.')[-1] == 'fits':
        save_votable = fitsfile.replace('.fits', '-sol.vot')
    else:
        save_votable = os.path.splitext(fitsfile)[-1] + '-sol.vot'    
        
    field = SkyCoord(head['RA']*u.deg, head['DEC']*u.deg)
    
    obs_end_time = Time(head['DATE-END']) #Time in 'isot'
    obs_time = Time(obs_end_time.iso, format='iso')
    try:
        head['TELESCOP']
    except(KeyError):
        head['TELESCOP'] = 'TAROT REUNION'
                
    if head['TELESCOP'] == 'TAROT CALERN':
        loc = '010'
        
    elif head['TELESCOP'] == 'TAROT CHILI':
        loc = '181'
        
    elif head['TELESCOP'] == 'TAROT REUNION':
        loc = '809'
    else:
        print("Not TAROT telescope")
        sys.exit()
        
    
    print("Observatory: %s\nObserved time: %s" %(loc, obs_time))
    print("Coordinate: %s" %field.to_string())
    
    # Search with the radius of 2k, TRE witll be crop to 2k, 4 images
    # Then search in each image like TCA and TCH.
    try:
        cobj = Skybot.cone_search(field,
                                  rad = 1.29*u.deg,
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
    
    if type(cobj) is QTable:
        cobj = Limit_Coor_to_Fits(fitsfile, catalog=cobj, ratext='RA', detext='DEC')
    else:
        pass
    
    if recode:
        if type(cobj) is QTable:
            cobj.write(save_votable, format='votable', overwrite=True)
    
    return cobj

def Check_Asteroid_Comet_Full_4k_FITS(fitsfile, recode=True):
        with fits.open(fitsfile) as hdu:
            head = hdu[0].header
            wcs = WCS(hdu[0].header)

        RA = head['RA']
        DEC = head['DEC']
        
        #       save solar object to votable
        if fitsfile.split('.')[-1] == 'new':
            save_votable = fitsfile.replace('.new', '-sol.vot')
        elif fitsfile.split('.')[-1] == 'fits':
            save_votable = fitsfile.replace('.fits', '-sol.vot')
        else:
            save_votable = os.path.splitext(fitsfile)[-1] + '-sol.vot'
            
        field = SkyCoord(RA*u.deg, DEC*u.deg)
        obs_end_time = Time(head['DATE-END']) #Time in 'isot'
        obs_time = Time(obs_end_time.iso, format='iso')
        
        try:
            head['TELESCOP']
        except(KeyError):
            head['TELESCOP'] = 'TAROT REUNION'
            
        if head['TELESCOP'] == 'TAROT CALERN':
            radius = 1.29
            loc = '010'
            print("%s %s" %(head['TELESCOP'], field.to_string()))
            cobj = Asteroid_Comet(field, radius, loc, obs_time)
            
        elif head['TELESCOP'] == 'TAROT CHILI':
            radius = 1.29
            loc = '181'
            print("%s %s" %(head['TELESCOP'], field.to_string()))
            cobj = Asteroid_Comet(field, radius, loc, obs_time)

        elif head['TELESCOP'] == 'TAROT REUNION':
            radius = 1.29
            loc = '809'

        #Divided to 4 part which equals to 2k but vstack in the end
            c1 = wcs.all_pix2world(1024, 1024, 0)
            field1 = SkyCoord(c1[0]*u.deg, c1[1]*u.deg)
            cobj1 = Asteroid_Comet(field1, radius, loc, obs_time)
            if type(cobj1) is QTable:
                pass
            else:
                cobj1 = ''
            
            c2 = wcs.all_pix2world(1024, 3072, 0)
            field2 = SkyCoord(c2[0]*u.deg, c2[1]*u.deg)
            cobj2 = Asteroid_Comet(field2, radius, loc, obs_time)
            if type(cobj2) is QTable:
                pass
            else:
                cobj2 = ''

            c3 = wcs.all_pix2world(3072, 3072, 0)
            field3 = SkyCoord(c3[0]*u.deg, c3[1]*u.deg)
            cobj3 = Asteroid_Comet(field3, radius, loc, obs_time)
            if type(cobj3) is QTable:
                pass
            else:
                cobj3 = ''
            
            c4 = wcs.all_pix2world(3072, 1024, 0)
            field4 = SkyCoord(c4[0]*u.deg, c4[1]*u.deg)
            cobj4 = Asteroid_Comet(field4, radius, loc, obs_time)
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
            
        print(head['TELESCOP'])
        print(obs_time)
        print("Coordinate: %s %s" %(RA, DEC))
        
        if type(cobj) is QTable:
            cobj = Limit_Coor_to_Fits(fitsfile, catalog=cobj, ratext='RA', detext='DEC')
        else:
            pass
    
        if recode:
            if type(cobj) is QTable:
                cobj.write(save_votable, format='votable', overwrite=True)
        return cobj
    
if __name__ == '__main__':
    home = os.path.expanduser('~')
    
    filefits2k = os.path.join(home, 'Documents/Data/Test/IM_20191207_023356000_000002_26943902.fits')
    filefits4k = os.path.join(home, 'Documents/Data/GrandMa_alert/IM_20190225_184237558_000001_76279408.fits')
    print(Check_Asteroid_Comet_Full_2k_FITS(filefits2k, 338.166967, -16.288820))
    
    print(Check_Asteroid_Comet_Full_4k_FITS(filefits4k))
    
    print(Check_Asteroid_Comet(filefits2k, 338.166967, -16.288820))
