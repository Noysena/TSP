#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 13:46:44 2019

@author: noysena
"""

import math
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from tarot.photometry.auto_photometry import Read_Count_Pixel

def Mag_Ref(fitsfile, star_flux):
    #Take a star at center image for reference star
    hdu = fits.open(fitsfile)
    header = hdu[0].header
    wcs = WCS(hdu[0].header)
    coor_candidate = SkyCoord(header['RA'], header['DEC'], frame = 'icrs',
                              unit = (u.deg, u.deg))
        
    catalog_retrieve = Vizier(column_filters={"Rmag":"< 16.0"},
                                  row_limit=-1).query_region(coor_candidate,
                                              radius=5*u.arcmin,
                                              verbose=False, catalog='NOMAD1')
        
    catalog_retrieve_nomad1 = catalog_retrieve[0]
        
    search_near_mag = catalog_retrieve_nomad1["Rmag"] < 15
    candidate_ref = catalog_retrieve_nomad1[search_near_mag]
    if len(candidate_ref) == 0:
        candidate_ref = catalog_retrieve_nomad1[0]["Rmag"]
    else:
        pass
    
    #Take the first star from table
    ref_RA = candidate_ref[0]["RAJ2000"]
    ref_DEC = candidate_ref[0]["DEJ2000"]
    ref_mag = candidate_ref[0]["Rmag"].item()
    
    Y, X = wcs.all_world2pix(ref_RA, ref_DEC, 0)
    Rflux, bgflx, bfstd, snr, stn = Read_Count_Pixel(fitsfile, X, Y)
    magnitude_ref = ref_mag + (-2.5*math.log10(star_flux/Rflux))
    return magnitude_ref

def Catalog_Nomad1(fitsfile):
    #Take a star at center image for reference star
    hdu = fits.open(fitsfile)
    header = hdu[0].header
    coor_candidate = SkyCoord(header['RA'], header['DEC'], frame = 'icrs',
                              unit = (u.deg, u.deg))
        
    catalog_retrieve = Vizier(column_filters={"Rmag":"< 17.0"},
                                  row_limit=-1).query_region(coor_candidate,
                                              radius=30*u.arcmin,
                                              verbose=False, catalog='NOMAD1')
        
    catalog_retrieve_nomad1 = catalog_retrieve[0]
        
    return catalog_retrieve_nomad1

def Mag_Ref_Table(fitsfile, tb_data):
    #data is in astropy table
    #Take a star at center image for reference star
    #update in table 'MAG_AUTO' and 'MAGERR_AUTO'
    #There is problem when the is no near reference star at RA & DEC
    hdu = fits.open(fitsfile)
    header = hdu[0].header
    wcs = WCS(hdu[0].header)
    
    coor_candidate = SkyCoord(header['RA'], header['DEC'], frame = 'icrs',
                              unit = (u.deg, u.deg))
        
    catalog_retrieve = Vizier(column_filters={"Rmag":"< 16.0"},
                                  row_limit=-1).query_region(coor_candidate,
                                              radius=5*u.arcmin,
                                              verbose=False, catalog='NOMAD1')
        
    catalog_retrieve_nomad1 = catalog_retrieve[0]
    search_near_mag = catalog_retrieve_nomad1["Rmag"] < 15
    candidate_ref = catalog_retrieve_nomad1[search_near_mag]
    #In case no ref star in clone search of 5*u.arcmin
    if len(candidate_ref) == 0:
        candidate_ref = catalog_retrieve_nomad1[0]["Rmag"]
    else:
        pass
    
    #Take the first star from table
    ref_RA = candidate_ref[0]["RAJ2000"]
    ref_DEC = candidate_ref[0]["DEJ2000"]
    ref_mag = candidate_ref[0]["Rmag"].item()
    
    Y, X = wcs.all_world2pix(ref_RA, ref_DEC, 0)
    Rflux, bgflx, bfstd, snr, stn = Read_Count_Pixel(fitsfile, X, Y)
    for i in range(len(tb_data)):
        star_mag = ref_mag + (-2.5*math.log10(tb_data[i]['FLUX_AUTO']/Rflux))
        erro_mag = ref_mag + (-2.5*math.log10(tb_data[i]['FLUXERR_AUTO']/Rflux))
        tb_data[i]['MAG_AUTO'] = np.round(star_mag, 4)
        tb_data[i]['MAGERR_AUTO'] = np.round(erro_mag, 4)
        
    return tb_data

if __name__ == "__main__":
    import os
    from astropy.table import Table
    home = os.path.expanduser('~')
    
    filefits = os.path.join(home, "Documents/Data/Test/IM_20191207_023356000_000002_26943902.fits")
    filecat = filefits.replace('.fits', '.cat')
    tbcat = Table.read(filecat, format='ascii.sextractor')
    
    mref = Mag_Ref(filefits, 20000)
    print('{:2.3f}'.format(mref))
    
    tbconvert = Mag_Ref_Table(filefits, tbcat)
    print(tbconvert)

    if False:
        with fits.open(filefits) as rf:
            head=rf[0].header
            data=rf[0].data
            wcs= WCS(rf[0].header)