#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 12:13:28 2019

@author: kanthanakorn
"""

#Search FITS in folder and check path

import os
import glob
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import Column
from sklearn.cluster import KMeans

home = os.path.expanduser('~')
Path_Test = os.path.join(home, "Documents/Data/GrandMa_alert")

# Brower path file
def List_File(Path):
    if os.access(Path, os.R_OK): #Ture
        findit = os.path.join(Path, '*.fits') 
        fitslist = sorted(glob.glob(findit))
        if not fitslist:
            findit = os.path.join(Path, '*.fits.gz')
            fitslist = sorted(glob.glob(findit))
    else:
        fitslist = None
        
    return fitslist

# Read FITS
def Fits_Info(File):
    hdu = fits.open(File)
    HDU_header = hdu[0].header
    wcs = WCS(File)
    hdu.close()
    return HDU_header, wcs

Image_list = List_File(Path_Test)
print(Image_list)

#%%
# Search FITS details DATE-OBS --> take only time
hdu_name = ["Triger-NAME", "FILENAME", "Time-OBS", "RA", "DEC", "COORDSYS",
            "FWHM", "EXPOSURE", "BGMEAN", "BGSIGMA", "FILTER", "TELESCOP"]

Max_FoV_TCA = 0
Max_FoV_TCH = 0
Max_FoV_TRE = 0
for i,j in enumerate(Image_list):
    hdu, wcs = Fits_Info(j)
    
    try:
        tele = hdu["TELESCOP"]
    except(KeyError):
        tele = "TAROT REUNION"
        
    if tele == "TAROT CALERN":
        try:
            FoV_TCA = int(hdu["NAME"][8:10])
        except(ValueError):
            try:
                FoV_TCA = int(hdu["NAME"][8:9])
            except(ValueError):
                FoV_TCA = 0
        else:
            pass
        
        if FoV_TCA > Max_FoV_TCA:
            Max_FoV_TCA = FoV_TCA
        
    if tele == "TAROT CHILI":
        try:
            FoV_TCH = int(hdu["NAME"][-1:])
        except(ValueError):
            try:
                FoV_TCH = int(hdu["NAME"][-2:])
            except(ValueError):
                FoV_TCH = 0
        else:
            pass
        
        if FoV_TCH > Max_FoV_TCH:
            Max_FoV_TCH = FoV_TCH
            
    if tele == "TAROT REUNION":
        try:
            FoV_TRE = int(hdu["NAME"][8:9])
        except(ValueError):
            try:
               FoV_TRE = int(hdu["NAME"][8:10])
            except(ValueError):
                FoV_TRE = 0
        else:
            pass
        
        if FoV_TRE > Max_FoV_TRE:
            Max_FoV_TRE = FoV_TRE
            
    try:
        fwhm = hdu["FWHM"]
    except(KeyError):
        fwhm = 99
            
    hdu_info = np.array([hdu["NAME"],
          hdu["FILENAME"],
          hdu["DATE-OBS"][11:23],
          np.around(hdu["RA"], decimals = 3),
          np.around(hdu["DEC"], decimals = 3),
          hdu["COORDSYS"],
          fwhm,
          hdu["EXPOSURE"],
          '{:0.3f}'.format(hdu["BGMEAN"]),
          '{:0.3f}'.format(hdu["BGSIGMA"]),
          hdu["FILTER"],
          tele])
#Image_table = Table(hdu_info, names=hdu_name)

    if i == 0:
        Image_table = Table(hdu_info, names=hdu_name,
                            masked=True,
                            meta={'name':'GW counterpart image'},
                            dtype=('S12', 'S45', 'S12', 'f8', 'f8', 'S4',
                                   'f8', 'f8', 'f8', 'f8', 'S1', 'S12'))
    else:
        Image_table.add_row(hdu_info)
        
Max_FoV = Max_FoV_TCA + Max_FoV_TCH + Max_FoV_TRE
   
Image_table["RA"].unit = u.deg
Image_table["DEC"].unit = u.deg
Image_table["FWHM"].unit = u.arcsec
Image_table["EXPOSURE"].unit = u.s

#print(Image_table)
print("Number of FoV : %d" %Max_FoV)
   
#%%
# Read FoV in FITS
def find_fov(fits):
    HDU_header, wcs = Fits_Info(fits)
    ra = round(HDU_header['CRVAL1'],4)
    dec = round(HDU_header['CRVAL2'],4)
    fov = [ra, dec]
    return fov

# Search each image for FoV
list_FOV=[]
for i in Image_list:
    FOV = find_fov(i)
    list_FOV.append(FOV)

FOV_of_data = np.array(list_FOV)    

if False:
    plt.figure()
    plt.plot([i[0] for i in FOV_of_data], [i[1] for i in FOV_of_data], 'g+')
    plt.title("FoV")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.grid(True)
    plt.show()
    
if True:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="mollweide")
    ax.plot([i[0] for i in FOV_of_data], [i[1] for i in FOV_of_data], 'g+')
    ax.grid(True)
    plt.show()
#%%
# Look for FoV in image's group
kmeans = KMeans(n_clusters=Max_FoV, init='k-means++', random_state=0).fit(FOV_of_data)

# Predict FoV and add to table
predicted_FoV = KMeans(n_clusters=Max_FoV, init='k-means++', random_state=0).fit_predict(FOV_of_data)
col_predicted_FoV = Column(predicted_FoV, name="FoV", dtype=int)
Image_table.add_column(col_predicted_FoV, index=0)

#Print Table of data
print(Image_table)


print ("TAROT Field of view\n %s" %kmeans.cluster_centers_)

if True:
    plt.figure()
    plt.plot([i[0] for i in FOV_of_data ], [i[1] for i in FOV_of_data], 'g+')
    plt.plot([j[0] for j in kmeans.cluster_centers_],[j[1] for j in kmeans.cluster_centers_], 'ro', markersize=10, markerfacecolor="None")
    plt.grid(True)
    plt.show()

    indir_path = os.path.dirname(Image_list[0])
    save_fov = os.path.join(indir_path, "FOVsData.npy")
    np.save(save_fov, kmeans.cluster_centers_)
    print ("List FOVs save in *.npy at : ", save_fov)

