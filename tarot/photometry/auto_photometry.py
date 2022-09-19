#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 12:24:59 2019

@author: kanthanakorn
"""

info = "Automated photometry after delete duplicate candidate in same FOV"""

import sys
import copy
import math
import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from tarot.utility.limit_in_fits import Limit_Coor_to_Fits

def Convert_W_P(ra, dec, hdu):
    wcs = WCS(hdu[0].header)
    x, y = wcs.all_world2pix(ra, dec, 0)
    return x, y

#Tiger time
def Triger_Date(trig_date, obs_date, time_format='hour'):
    #event format isof : 2019-04-12T05:30:44.165622
    if isinstance(trig_date, str):
        pass
    else:
        sys.exit(0)
        
	 #Triger by LIGO or the first observation with TAROT
    trig_dpt = pd.to_datetime(trig_date.replace('T', ' '))
    
    #Observation time
    obs_dtp = pd.to_datetime(obs_date.replace('T' , ' '))
    
    delta_dtp = obs_dtp - trig_dpt
    #Provide time in Timedelta after triger
    if time_format == 'second':
#    delta_dtp = ("%s" %(obs_dtp - trig_dpt))
        delta_dpt_out = delta_dtp.total_seconds()
    elif time_format == 'minute':
        delta_dpt_out = delta_dtp.total_seconds() / 360.
    elif time_format == 'hour':
        delta_dpt_out = delta_dtp.total_seconds() / 3600.
    elif time_format == 'day':
        delta_dpt_out = delta_dtp.total_seconds() / (3600.*24)
    else:
        delta_dpt_out = delta_dtp.total_seconds()
    return delta_dpt_out

def Next_Circle(sep_ref_list, limit_angle):
    #limit angle in arcsecond
    if len(sep_ref_list) < 1:
        candi_ref = None
    else:
        for i in sep_ref_list:
            if i > limit_angle:
                candi_ref = i
                break
            else:
                candi_ref = None
 
    return candi_ref

def Retrieve_Ref(RA, DEC): #, CenterRA, CenterDEC):
        
    coor_candidate_fk5 = SkyCoord(RA, DEC, frame = 'icrs',
                                           unit = (u.deg, u.deg))
        
    catalog_retrieve = Vizier(column_filters={"Rmag":"< 17.0"},
                                  row_limit=-1).query_region(coor_candidate_fk5,
                                              radius=10*u.arcmin,
                                              verbose=False, catalog='NOMAD1')
    
    catalog_retrieve_nomad1 = catalog_retrieve[0]
#    print(catalog_retrieve_nomad1)
    search_near_mag = catalog_retrieve_nomad1["Rmag"] < 15
    if np.any(search_near_mag):
        candidate_ref = catalog_retrieve_nomad1[search_near_mag]
    else:
        search_near_mag = catalog_retrieve_nomad1["Rmag"] < 16
        candidate_ref = catalog_retrieve_nomad1[search_near_mag]
        if np.any(search_near_mag):
            candidate_ref = catalog_retrieve_nomad1
        
    
    coor_candidate_ref_fk5 = SkyCoord(candidate_ref["RAJ2000"],
                                  candidate_ref["DEJ2000"],
                                  frame = 'fk5', unit = (u.deg, u.deg))
    
    coor_candidate_ref = coor_candidate_ref_fk5.transform_to('icrs')
    
    sep_candi_ref = coor_candidate_fk5.separation(coor_candidate_ref)
    
    sort_ref = sorted(sep_candi_ref.arcsec)
    
    #There was the case that no ref near by, Need other methods
    min_sep_ref = Next_Circle(sort_ref, 80.0)
    if min_sep_ref == None:
        min_sep_ref = Next_Circle(sort_ref, 60.0)
    else:
        pass
    
    sellect_min_sep_ref = np.where(min_sep_ref == np.array(sep_candi_ref.arcsec))
    Nomad1_reference = candidate_ref[sellect_min_sep_ref]
    
    # Reference checking
    if False:
        print(sep_candi_ref.arcsec)
        print(sort_ref)
        print("Choose separation Ref[arcsec]: ", min_sep_ref)
        print(sellect_min_sep_ref)
        print(candidate_ref[sellect_min_sep_ref])
        print("Ref coor: ",
              Nomad1_reference["RAJ2000"],
              Nomad1_reference["DEJ2000"],
              Nomad1_reference["Rmag"])
        
        min_coor_candidate_ref_fk5 = SkyCoord(Nomad1_reference["RAJ2000"],
                                          Nomad1_reference["DEJ2000"],
                                          frame = 'fk5',
                                          unit = (u.deg, u.deg))
        
        min_coor_candidate_ref = min_coor_candidate_ref_fk5.transform_to('icrs')

        check_coor = coor_candidate_fk5.separation(min_coor_candidate_ref)
        
        print("Candi coor: ", RA, DEC)
        print("In table; I take [arcsec]: ", check_coor.arcsec)
        print("This is the closest Ref: ", check_coor.arcsec == min_sep_ref)
    
    return Nomad1_reference

def Retrieve_Ref_2(filefits, RA, DEC):
        
    coor_candidate_fk5 = SkyCoord(RA, DEC, frame = 'icrs',
                                           unit = (u.deg, u.deg))
        
    catalog_retrieve = Vizier(column_filters={"Rmag":"< 17.0"},
                                  row_limit=-1).query_region(coor_candidate_fk5,
                                              radius=10*u.arcmin,
                                              verbose=False, catalog='NOMAD1')
    #Limit catalog in FITS
    catalog_retrieve_nomad1 = Limit_Coor_to_Fits(filefits=filefits, catalog=catalog_retrieve[0])

    search_near_mag = catalog_retrieve_nomad1["Rmag"] < 15
    if np.any(search_near_mag):
        candidate_ref = catalog_retrieve_nomad1[search_near_mag]
    else:
        search_near_mag = catalog_retrieve_nomad1["Rmag"] < 16
        candidate_ref = catalog_retrieve_nomad1[search_near_mag]
        if np.any(search_near_mag):
            candidate_ref = catalog_retrieve_nomad1
        
    
    coor_candidate_ref_fk5 = SkyCoord(candidate_ref["RAJ2000"],
                                  candidate_ref["DEJ2000"],
                                  frame = 'fk5', unit = (u.deg, u.deg))
    
    coor_candidate_ref = coor_candidate_ref_fk5.transform_to('icrs')
    
    sep_candi_ref = coor_candidate_fk5.separation(coor_candidate_ref)
    
    sort_ref = sorted(sep_candi_ref.arcsec)
    
    #There was the case that no ref near by, Need other methods
    min_sep_ref = Next_Circle(sort_ref, 80.0)
    if min_sep_ref == None:
        min_sep_ref = Next_Circle(sort_ref, 60.0)
    else:
        print("No reference source is next to target: %s %s." %(RA, DEC))
    
    sellect_min_sep_ref = np.where(min_sep_ref == np.array(sep_candi_ref.arcsec))
    Nomad1_reference = candidate_ref[sellect_min_sep_ref]
    
    # Reference checking
    if False:
        print(sep_candi_ref.arcsec)
        print(sort_ref)
        print("Choose separation Ref[arcsec]: ", min_sep_ref)
        print(sellect_min_sep_ref)
        print(candidate_ref[sellect_min_sep_ref])
        print("Ref coor: ",
              Nomad1_reference["RAJ2000"],
              Nomad1_reference["DEJ2000"],
              Nomad1_reference["Rmag"])
        
        min_coor_candidate_ref_fk5 = SkyCoord(Nomad1_reference["RAJ2000"],
                                          Nomad1_reference["DEJ2000"],
                                          frame = 'fk5',
                                          unit = (u.deg, u.deg))
        
        min_coor_candidate_ref = min_coor_candidate_ref_fk5.transform_to('icrs')

        check_coor = coor_candidate_fk5.separation(min_coor_candidate_ref)
        
        print("Candi coor: ", RA, DEC)
        print("In table; I take [arcsec]: ", check_coor.arcsec)
        print("This is the closest Ref: ", check_coor.arcsec == min_sep_ref)
    
    return Nomad1_reference

def Read_Count_Pixel(fitsfile, Xpoint, Ypoint):
    #Read fits
    hdu = fits.open(fitsfile)
    flux = hdu[0].data
    hdu.close()
    # state : source area, background area, CCD gain
    Gain = 1.0
    source_area = 10
    bg_area = 20
    nbpix = ((source_area*2))**2

    #Source
    Xmin = int(Xpoint - source_area)
    Ymin = int(Ypoint - source_area)

    Xmax = int(Xpoint + source_area)
    Ymax = int(Ypoint + source_area)

    #Copy area that contains source
    source = copy.copy(flux[Xmin:Xmax, Ymin:Ymax])

    #Background
    BGXmin = int(Xpoint - bg_area)
    BGYmin = int(Ypoint - bg_area)

    BGXmax = int(Xpoint + bg_area)
    BGYmax = int(Ypoint + bg_area)
    
    BGXmin2 = bg_area - source_area
    BGYmin2 = bg_area - source_area
    
    BGXmax2 = bg_area + source_area
    BGYmax2 = bg_area + source_area

    #Find background median from full image
    BG_ = copy.copy(flux[BGXmin:BGXmax, BGYmin:BGYmax])

    for counter in range(3):
        x=0
        while x <= 20:
    #Obtain median, std and peak to eliminate from data
            BG_stat = np.median(BG_)
            BG_std = np.round(BG_.std(), 5)
            Upper_limit = BG_stat + BG_.std() * 5
    
            cut_arr = np.where(BG_ > Upper_limit)
            BG_[cut_arr] = Upper_limit
        
            BG_std_after_limit = np.round(BG_.std(), 5)
            
            
            if BG_std_after_limit == BG_std:
                break
            x += 1
    
    #Take out source from counting area and replace with BG_stat
        BG_[BGXmin2:BGXmax2,BGYmin2:BGYmax2] = float(BG_stat)
    
    BG_stat =  np.median(BG_)*Gain
    BG_std = BG_.std()*Gain    
    
#    source = copy.copy(flux[Xmin:Xmax, Ymin:Ymax])
    source_matrix = (source) * Gain
    
    source_stat = np.sum(source_matrix) - nbpix*BG_stat
    source_max = np.max(source_matrix)
    source_signal = source_max - BG_stat
    noise = BG_std
    snr = source_signal / noise
    
    #If Signol to Nose is less than 2 than there is no source in image
    if snr < 2.:
        source_stat = 1.

    if False:
        print("Candidate ##################################### Reference flux")
        print('Sum of source : {: 0.2f}'.format(np.sum(source_matrix)))
        print('nbpix : {:d}'.format(nbpix))
        print('Flux : {: 0.2f}'.format(source_stat))
        print('Background : {: 0.2f}'.format(BG_stat))
        print('Background STD : {: 0.2f}'.format(BG_std))
        print('SNR : {: 0.2f}'.format(snr))
        
    if False:
        fluxmin = np.median(source)-500
        fluxmax = np.median(source)+500
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(131, aspect='equal')
        ax1.imshow(source,
                   vmin=fluxmin, vmax=fluxmax, origin='lower', cmap='Greys')
    
        ax2 = fig1.add_subplot(132, aspect='equal')
        ax2.imshow(flux,
                   vmin=fluxmin, vmax=fluxmax, origin='lower', cmap='Greys')
        ax2.plot(Ypoint, Xpoint, 'r+')
    
        ax2.plot(Ypoint+30, Xpoint, 'b+')
        ax2.plot(Ypoint, Xpoint + 30, 'g+')
    
        ax2.plot(Ypoint-30, Xpoint, 'b+')
        ax2.plot(Ypoint, Xpoint - 30, 'g+')

        ax3 = fig1.add_subplot(133, aspect='equal')
        ax3.imshow(BG_, 
                   vmin=fluxmin, vmax=fluxmax, origin='lower', cmap='Greys')
    
    return source_stat, BG_stat, BG_std, snr, source_signal

def Run_Photrometry(fitsfile, ra, dec):
    #Read fits
    hdu = fits.open(fitsfile)
    header = hdu[0].header
    observe_date = header["DATE-OBS"]
    wcs = WCS(fitsfile)
    hdu.close()
    
    # Call Reference star in case we have only one
    ref = Retrieve_Ref(fitsfile, ra, dec)
    ref_ra = ref["RAJ2000"]
    ref_dec = ref["DEJ2000"]
    ref_mag = ref["Rmag"]
    
    #Obtain pixel coordinate from WCS
    Y, X = wcs.all_world2pix(ra, dec, 0)
    ref_Y, ref_X = wcs.all_world2pix(ref_ra, ref_dec, 0)
    
    #Phototmetry candidate if flux is negative, then candidate flux is 1
    flux_candi, flux_candi_BG, candi_STD, candi_SNR, flux_central_candi = Read_Count_Pixel(fitsfile, X, Y)
    if flux_candi < 0:
        flux_candi = 1
    
    #Photmetry reference if flux is negative, then candidate flux is 1
    flux_refer, flux_refer_BG, refer_STD, refer_SNR, flux_central_refer = Read_Count_Pixel(fitsfile, ref_X, ref_Y)
    if flux_refer < 0:
        flux_refer = 1

    # If Signal to Nosie is more than 2 then normally computation
    if flux_candi > 1 and flux_refer > 1 and candi_SNR > 2. :
        candi_mag = ref_mag + (-2.5*math.log10(flux_candi/flux_refer))
    
    #If Signal to Noise is less than 2 then candidate flux is 3STD of Back ground reference star
    rho = flux_refer / flux_central_refer
    limiting_mag = ref_mag + (-2.5*math.log10(3*refer_STD*rho/flux_refer))
    if flux_candi == 1 and flux_refer > 1 :
        candi_mag = limiting_mag
        
    #This case is the reference star is bad but flux_candi is measuarable.
    if flux_candi >= 1 and flux_refer == 1:
        candi_mag = limiting_mag
        
    raw_ref_mag = -2.5*math.log10(flux_refer)
    offset_mag = ref_mag - raw_ref_mag
    candi_deltamag = 1. / candi_SNR
    
    if flux_candi == 'nan':
        datalist = 0
    else:
        datalist = [np.round(flux_candi,3), np.round(flux_candi_BG, 3),
                np.round(candi_STD, 3), np.round(candi_SNR, 3),
                np.round(flux_refer,3), np.round(flux_refer_BG, 3),
                np.round(refer_STD, 3), np.round(refer_SNR, 3),
                np.round(offset_mag,3), np.round(limiting_mag, 3),
                np.round(candi_mag, 3), np.round(candi_deltamag, 3),
                observe_date]

    return datalist


def Run_Photrometries(fitsfile, ra, dec, ref_ra, ref_dec, ref_mag):
    #Read fits
    hdu = fits.open(fitsfile)
    header = hdu[0].header
    observe_date = header["DATE-OBS"]
    wcs = WCS(fitsfile)
    hdu.close()
    
    #Obtain pixel coordinate from WCS
    Y, X = wcs.all_world2pix(ra, dec, 0)
    ref_Y, ref_X = wcs.all_world2pix(ref_ra, ref_dec, 0)
    
    #Phototmetry candidate if flux is negative, then candidate flux is 1
    try:
        flux_candi, flux_candi_BG, candi_STD, candi_SNR, flux_central_candi = Read_Count_Pixel(fitsfile, X, Y)
        A = True
        if flux_candi < 0:
            flux_candi = 1
        else:
            pass
    except:
        A = False
    
    #Photmetry reference if flux is negative, then candidate flux is 1
    try:
        flux_refer, flux_refer_BG, refer_STD, refer_SNR, flux_central_refer = Read_Count_Pixel(fitsfile, ref_X, ref_Y)
        B = True
        if flux_refer < 0:
            flux_refer = 1
        else:
            pass
    except:
        B = False

    if A and B:
        # If Signal to Nosie is more than 2 then normally computation
        if flux_candi > 1 and flux_refer > 1 and candi_SNR > 2. :
            candi_mag = ref_mag + (-2.5*math.log10(flux_candi/flux_refer))
        else:
            pass
            #If Signal to Noise is less than 2 then candidate flux is 3STD of Back ground reference star
        rho = flux_refer / flux_central_refer
        limiting_mag = ref_mag + (-2.5*math.log10(3*refer_STD*rho/flux_refer))
        
        if flux_candi == 1 and flux_refer > 1 :
            candi_mag = limiting_mag
        else:
            pass
        
        #This case is the reference star is bad but flux_candi is measuarable.
        if flux_candi >= 1 and flux_refer == 1:
            candi_mag = limiting_mag
        else:
            pass
        
        raw_ref_mag = -2.5*math.log10(flux_refer)
        offset_mag = ref_mag - raw_ref_mag
        candi_deltamag = 1. / candi_SNR
    
        if flux_candi == 'nan':
            datalist = [None, None, None, None, None, None, None,
                        None, None, None, None, None, None]
        else:
            datalist = [np.round(flux_candi,3), np.round(flux_candi_BG, 3),
                        np.round(candi_STD, 3), np.round(candi_SNR, 3),
                        np.round(flux_refer,3), np.round(flux_refer_BG, 3),
                        np.round(refer_STD, 3), np.round(refer_SNR, 3),
                        np.round(offset_mag,3), np.round(limiting_mag, 3),
                        np.round(candi_mag, 3), np.round(candi_deltamag, 3),
                        observe_date]
    else:
        datalist = [None, None, None, None, None, None, None,
                        None, None, None, None, None, None]
    return datalist

if __name__ == '__main__':
    from tarot.photometry.star_ref import Catalog_Nomad1

    filefits = '/home/tarot/Documents/Data/GrandMa_alert/IM_20190225_184237558_000001_76279408.fits'
    with fits.open(filefits) as rf:
        hdu = rf[0].header
        
    RA = hdu['RA']
    DEC = hdu['DEC']

    Nomad1_reference = Retrieve_Ref_2(filefits, RA, DEC)
    print(Nomad1_reference)

    Nomad1_ref = Catalog_Nomad1(filefits)
    print(Nomad1_ref)
