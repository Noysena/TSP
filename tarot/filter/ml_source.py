#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:19:33 2019

@author: noysena

ML to identify sources
"""

import numpy as np
from joblib import load
from astropy.wcs import WCS
from astropy.io import fits


class ML_SOURCE(object):
    
    def MAIN():
        return("Machine learning to identify source.")
        
    def __init__(self, file_sif):
        self.clf = load(file_sif)
        self.predicted = ''
        
    def What_Source(self, image):
        image = image.reshape(1,400)
        self.predicted = self.clf.predict(image)
        
if __name__ == '__main__':
    import os
    from skimage import exposure
    import matplotlib.pyplot as plt
    from tarot.utility.general_funtion import sub_image
    from tarot.utility.general_funtion import complet_array
    from astropy.table import Table
    
    run = ML_SOURCE('source_identification_4k.joblib')
    
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
            print("%d : %s" %(i, img_contrast.shape))
            run.What_Source(source)
#            print("source: %d" %run.predicted)
            
            plt.figure()
            plt.imshow(source, cmap = 'gray')
            plt.title("Source %s" %run.predicted)
            plt.show()
        
            run.What_Source(source)
            print("source: %d" %run.predicted)
    else:
        pass
            
    if False:
        home = os.path.expanduser('~')
        filedat = os.path.join(home, 'Documents/Data/Test/IM_20191207_023356000_000002_26943902.candi.dat')
        filefits = filedat.replace('.candi.dat', '.fits')
        
        with fits.open(filefits) as rf:
            hdu = rf[0].header
            wcs = WCS(rf[0].header)
            data = rf[0].data
            
        tbdat = Table.read(filedat, format='ascii')
        print(tbdat)
        keep = []
        for i, j in enumerate(tbdat):
            source = sub_image(j['Y_IMAGE'], j['X_IMAGE'],  10, data)
            p2, p98 = np.percentile(source, (15, 85))
            img_contrast = exposure.rescale_intensity(source, in_range=(p2, p98))
            img_contrast = complet_array(img_contrast, 20)
            
            plt.figure()
            plt.imshow(source, cmap = 'gray')
            plt.show()
        
            run.What_Source(source)
            print("source: %d" %run.predicted)
            if run.predicted == 3:
                keep.append(i)
                
            else:
                pass
            
            if run.predicted == 0:
                keep.append(i)
            else:
                pass
           
        tbdat.remove_rows(keep)
        print(tbdat)
    else:
        pass
        