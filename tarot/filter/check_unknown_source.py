#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:21:21 2018

@author: kanthanakorn
"""
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
# =============================================================================
# Search unknow source in USNO-B1.o catalog in KRON_RADIUS   
# I use the X_WORLD & Y_WORLD convert form pixel + SIP5
# RA and DEC form iteration result will be ignored as we want to refer back to
# source by read form FITS WCS
# =============================================================================
def Check_Object(Data_candi, Cdata_candi):
    candidate_idx = []
    for i in range(len(Cdata_candi)):
        #Build radius for flexible search in USNO-B1.0
        Kron_radius = float(Data_candi["KRON_RADIUS"][i])
        Search_radius = Kron_radius * 4.0
        
        """ See what happen in checking step """
        Want_to_see = False
        if Want_to_see:
            print("Search radius in USNO-B1.0(arcsec): ", str(Search_radius))
        else:
            pass
            
#        confirm_candi0 = LoadCata_HIP(Cdata[i].ra, Cdata[i].dec,16.0)
        confirm_candi1 = Vizier(column_filters={"R1mag":"< 18"},
                                row_limit=-1).query_region(Cdata_candi[i],
                                            radius=Search_radius*u.arcsec,
                                            catalog = "I/284",
                                            verbose=False,
                                            cache=False)
        
        print("Searching for candidates..." ,end = '\r')
        
        if str(confirm_candi1) != "Empty TableList":
            USNO_in_radius = SkyCoord(confirm_candi1[0]["RAJ2000"],
                                      confirm_candi1[0]["DEJ2000"],
                                      frame = 'icrs',
                                      unit = (u.deg, u.deg))
            AngSept = USNO_in_radius.separation(Cdata_candi[i])
        else:
            pass
           
            if Want_to_see:
                print(AngSept.arcsec)
                print(confirm_candi1[0]["R1mag"])
            else:
                pass
            
        if str(confirm_candi1) == "Empty TableList":
            candidate_idx.append(i)
        else:
            pass
    return candidate_idx

def Check_Object_NOMAD(Data_candi, Cdata_candi):
    candidate_idx = []
    for i in range(len(Cdata_candi)):
        #Build radius for flexible search in USNO-B1.0
        Kron_radius = float(Data_candi["KRON_RADIUS"][i])
        Search_radius = Kron_radius * 4.0
        
        """ See what happen in checking step """
        Want_to_see = False
        if Want_to_see:
            print("Search radius in USNO-B1.0(arcsec): ", str(Search_radius))
        else:
            pass
            
#        confirm_candi0 = LoadCata_HIP(Cdata[i].ra, Cdata[i].dec,16.0)
        confirm_candi1 = Vizier(column_filters={"Rmag":"< 20"},
                                row_limit=-1).query_region(Cdata_candi[i],
                                            radius=Search_radius*u.arcsec,
                                            catalog = "I/297",
                                            verbose=False,
                                            cache=False)
        
        print("Searching for candidates..." ,end = '\r')
        
        if str(confirm_candi1) != "Empty TableList":
            USNO_in_radius = SkyCoord(confirm_candi1[0]["RAJ2000"],
                                      confirm_candi1[0]["DEJ2000"],
                                      frame = 'icrs',
                                      unit = (u.deg, u.deg))
            AngSept = USNO_in_radius.separation(Cdata_candi[i])
        else:
            pass
           
            if Want_to_see:
                print(AngSept.arcsec)
                print(confirm_candi1[0]["R1mag"])
            else:
                pass
            
        if str(confirm_candi1) == "Empty TableList":
            candidate_idx.append(i)
        else:
            pass
    return candidate_idx
