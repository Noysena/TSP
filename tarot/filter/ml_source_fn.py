#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 17:46:37 2020

@author: kanthanakorn
"""
import os
from joblib import load

def Id_Source(image, header):
    if header['TELESCOP'] == 'TAROT CALERN':
        file_sif = 'source_identification_2k.joblib'
    elif header['TELESCOP'] == 'TAROT CHILI':
        file_sif = 'source_identification_2k.joblib'
    elif header['TELESCOP'] == 'TAROT REUNION':
        file_sif = 'source_identification_4k.joblib'
    else:
        print("Not TAROT telescope")
    home = os.path.expanduser('~')
    sif = os.path.join(home, 'TSP/tarot/filter/' + file_sif)
    clf = load(sif)
    ima = image.reshape(1,400)
    predicted = clf.predict(ima)
    
    if predicted == 0:
        source = False
    elif predicted == 1:
        source = True
    elif predicted == 2:
        source = True
    elif predicted == 3:
        source = False
    elif predicted == 4:
        source = False
    else:
        source = False
    
    return source

if __name__ == '__main__':
    import numpy as np
    from skimage import exposure
    from astropy.wcs import WCS
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from tarot.utility.general_funtion import sub_image
    from tarot.utility.general_funtion import complet_array
#    from astropy.table import Table
    
    
    
    if True:
        #Need to build predicted for TRE image.
        home = os.path.expanduser('~')
        filefits = os.path.join(home, 'Documents/Data/GrandMa_alert/IM_20190225_184705558_000001_76279408.fits')
        
        with fits.open(filefits) as rf:
            hdu = rf[0].header
            wcs = WCS(rf[0].header)
            data = rf[0].data
        
        #import data in x,y pixel, and it swap between  x and y
        sl = {1:[2048, 2016], 2:[2062, 2034], 3:[2071, 2039], 4:[2047, 2056],
              5:[2014, 2064], 6:[2071, 2048], 7:[2031, 2073], 8:[1995, 2074],
              9:[3758, 1616], 10:[3729, 1634], 11:[3831, 1638], 12:[3810, 1586]}
        for i, j in sl.items():
            source = sub_image(j[1], j[0], 10, data)
            p2, p98 = np.percentile(source, (15, 85))
            img_contrast = exposure.rescale_intensity(source, in_range=(p2, p98))
    #        img_contrast = exposure.equalize_adapthist(source)
            img_contrast = complet_array(img_contrast, 20)
            
            run = Id_Source(img_contrast, hdu)
            print("%d : %s" %(i, img_contrast.shape))
            
            plt.figure()
            plt.imshow(source, cmap = 'gray')
            plt.title("Source: %s" %run)
            plt.show()
        
            
    else:
        pass