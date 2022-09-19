#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 11:41:17 2020

@author: kanthanakorn
"""
import os
import time
import subprocess

def Astrometry(ra, dec, dir_fitsfile, fitsfile):
    fitsnew = fitsfile.replace('.fits', '.new')
    Radius = 6.0
    a=0
    runcode = ("solve-field --fits-image --guess-scale --no-plots\
                       --skip-solved --corr none --index-xyls none --match none\
                       --solved none --rdls none --overwrite --tweak-order 5 \
                       --ra %s --dec %s --radius %s --dir %s %s" 
                       %(ra, dec, Radius, dir_fitsfile, fitsfile))
        
    subprocess.Popen(runcode, stdout=subprocess.PIPE, shell=True)
    
    while not os.path.exists(fitsnew):
            a += 1
            print("Field is solving...%s" %(a), end='\r')
            time.sleep(1)
            
            if a == 200:
                print("Bad image, No callibration result.")
                break