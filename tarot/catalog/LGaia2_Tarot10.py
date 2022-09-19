#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Edited : K Noysena, noysena dot kanthanakorn dot oca dot eu
# PhD student: Aritemis
# Original code: Thomas Boch (github: tboch) - thomas dot boch at astro dot unistra dot fr
# Affiliation: Observatoire astronomique de Strasbourg


"""
This will download data directly from achive in Tarot 10
"""
import sys, os
import csv
import json
import math
import healpy
from astropy.table import Table
#import timeit

#Read error if there is exist.
def output_error(msg, exit=True):
    print('Error name: %s' %msg)
    if exit:
        sys.exit()

#build radius for input coordinate   
def radec2thetaphi(ra, dec):
    return math.radians((90-dec)), math.radians(ra) 

#define location for database; here is 'cgi-config.json'
def get_cgi_config_file_name():
    return 'cgi-config.json'

#obtain database of Gaia from folder; pointing by FN 'get_cgi_config_file_name'
def get_metafile_path(root):
    return os.path.join(root, 'metadata.json')

def get_path(root, nside, ipix):
    """
    Return path of file for given root directory,
    nside and ipix
    """
    dir_idx = (ipix // 10000)*10000;
    return os.path.join(root, "nside%d/dir%d/npix%d.csv" % (nside, dir_idx, ipix))    

def fields_as_votable(fields):
    sb = []
    for f in fields:
        nu = ("%s," %f['name'])
        sb.append(nu) #Take .encolde out to down form class list to type list
    return sb +['\n']

def sph_dist(ra1, dec1,ra2, dec2):
    """
    Compute the spherical distance between 2 pairs of coordinates
    using the Haversine formula
    
    Input coordinates are in decimal degrees
    Output: angular distance in decimal degrees
    """
    ra1_rad  = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad  = math.radians(ra2)
    dec2_rad = math.radians(dec2)

    d = math.sin((dec1_rad-dec2_rad)/2)**2;
    d += math.sin((ra1_rad-ra2_rad)/2)**2 * math.cos(dec1_rad)*math.cos(dec2_rad)

    return math.degrees(2*math.asin(math.sqrt(d)))

def LoadCata(RA, DEC, SR):
    #RA, DEC and SR in deg
    #check if config file is present and use the same as cs.py to make sure that *.json file is there
#    script_dir = os.path.abspath(os.path.dirname("/mnt/ext_disk/cgi-bin/cs.py"))
    script_dir = os.path.abspath(os.path.dirname("/mnt/tarot/gaia_data/cgi-bin/cs.py"))
    conf_path = os.path.join(script_dir, get_cgi_config_file_name())
    data_path = None

    if not os.path.exists(conf_path):
    # perhaps data is in the same directory
        metadata_path = get_metafile_path('.')
        if not os.path.exists(metadata_path):
            output_error('Service error: could not find config file %s' % (conf_path))
        else:
            data_path = os.path.abspath('.')
    else:
        with open(conf_path) as h:
            config = json.loads(h.read())
    
        data_path = os.path.abspath(config['dataPath'])
    
    metadata_path = get_metafile_path(data_path)
    if not os.path.exists(metadata_path):
        output_error('Service error: could not find metadata file %s' % (metadata_path))
    
    with open(metadata_path) as h:
        metadata = json.loads(h.read())
        
# retrieve info
    fields = metadata['fields']
    nside = metadata['nside']


#def loadCata(RA,DEC,SR):
    try:
        ra = float(RA)
    except:
        output_error("Could not parse value '%f' of RA parameter as a float" % (RA))
    try:
        dec = float(DEC)
    except:
        output_error("Could not parse value '%f' of DEC parameter as a float" % (DEC))
    try:
        sr = float(SR)
    except:
        output_error("Could not parse value '%f' of SR parameter as a float" % (SR))

# Check if parameters are withing sensible range
    if ra<0 or ra>=360:
        output_error('Value for RA parameter should be in range [0, 360[')
    if dec<-90 or dec>90:
        output_error('Value for DEC parameter should be in range [-90, 90]')
    if sr<0:
        output_error('value for SR parameter should be >=0')

    ra_idx = None
    dec_idx = None
    k = -1
    for f in fields:
        k += 1
        if ra_idx and dec_idx:
            break
        if 'ucd' in f:
            if f['ucd']=='POS_EQ_RA_MAIN':
                ra_idx = k
            elif f['ucd']=='POS_EQ_DEC_MAIN':
                dec_idx = k
    
    if ra_idx==None:
        output_error("Could not find field with ucd='POS_EQ_RA_MAIN'. Missing info in %s" % (metadata_path))
    if dec_idx==None:
        output_error("Could not find field with ucd='POS_EQ_DEC_MAIN'. Missing info in %s" % (metadata_path))
        
#where to save file, we will put in /tmp
    cata = open('/tmp/cata_GAIA.txt', 'w')
    data_name = fields_as_votable(fields)
    cata.writelines("%s" %nom for nom in data_name)
#    print(type(data_name), data_name)

# healpix query to retrieve data in requested cone
    keep = []
    theta, phi = radec2thetaphi(ra, dec) 
    vec = healpy.ang2vec(theta, phi)
    healpix_cells = healpy.query_disc(nside, vec, math.radians(sr), inclusive=True, nest=True)
#    print(healpix_cells)

    for ipix in healpix_cells:
        ipix_path = get_path(data_path, nside, ipix)

        if not os.path.exists(ipix_path):
            continue
        with open(ipix_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
        # read data in csv and write to ascii
            for row in reader:        
                row_ra  = float(row[ra_idx])
                row_dec = float(row[dec_idx])
                dist = sph_dist(ra, dec, row_ra, row_dec)
                if dist>sr:
                    continue
                add_line=','.join(row)+'\n'
                keep.append(add_line)
#    print(keep)
    cata.writelines("%s" %donner for donner in keep)
    cata.close()

#Need to recode without writing it down to hardisk
#read ascii and then convert to Table from, not request by TSP    
    Read_ascii_cata = open('/tmp/cata_GAIA.txt', 'r')
    Read_cata = Read_ascii_cata.read()
    Read_ascii_cata.close()

#Convert  ASCII read to Table from CSV, not request by TSP
    LCatalog  = Table.read(Read_cata, format='csv')
    Read_ascii_cata.close()

    return LCatalog

if __name__ == '__main__':
    home = os.path.expanduser('~')
    fitsfile = os.path.join(home, "Documents/Data/Test/IM_20191207_023356000_000002_26943902.fits")
    from astropy.io import fits
    with fits.open(fitsfile) as rfit:
        hd = rfit[0].header
        
    LGaia = LoadCata(hd['RA'], hd['DEC'], 2.0)
    LGaia.write(fitsfile.replace(".fits", ".cata.csv"), overwrite=True)
