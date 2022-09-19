#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 22:34:00 2019

@author: noysena
"""

import os
import sys
import glob
import json
import numpy as np
from astropy.io import fits

info = "Read FITS and export dictionary and dump to json format"

folder = os.path.join(sys.argv[1], "*.fits")
#folder = os.path.join(sys.argv[1], "*.new")
filelist = glob.glob(folder)
gws={}
for i in filelist:
    with fits.open(i) as hdul:
        hd = hdul[0].header
        center = {hd["FILENAME"] : {"FOV" : 999,
				    "RA" : np.round(hd["RA"], 4),
				    "DEC" : np.round(hd["DEC"], 4),
				    "SNAME": hd["SNAME"],
				    "DATE-OBS" : hd["DATE-OBS"]}}
    print(center)
    gws.update(center)

#Check FOV and rewrite json and save FOV
save_name = os.path.basename(sys.argv[1]) + '_FOV.txt'
textname = os.path.join(sys.argv[1], save_name)
mylist = []
for i, j in gws.items():
	list_txt = ("%s_%s" %(gws[i].get('RA'), gws[i].get('DEC')))
	mylist.append(list_txt)

#Eliminate duplicated FOV, then we have tiles of observation.
mylist = list(dict.fromkeys(mylist))
print("Tiles No. : %s" %len(mylist))

FOV = {}
for i, j in enumerate(mylist):
	ra, dec = j.rsplit("_")
	FOV.update({i: {'RA' : ra, 'DEC' : dec}})
print(FOV)
dumpfov = json.dumps(FOV)
with open(textname, 'w') as fov:
	fov.write(dumpfov)

#Update FOV in dict
for i, j in FOV.items():
	for m, n in gws.items():
		if float(gws[m].get('RA')) == float(FOV[i].get('RA')) and float(gws[m].get('DEC')) == float(FOV[i].get('DEC')):
			gws[m].update({'FOV' : i})
		else:
			pass
#Build JSON and write to file
dumpjs = json.dumps(gws)
save_name = os.path.basename(sys.argv[1]) + '.txt'
textname = os.path.join(sys.argv[1], save_name)
with open(textname, 'w') as rj:
	rj.write(dumpjs)