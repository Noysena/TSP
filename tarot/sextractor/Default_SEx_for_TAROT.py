#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Created on Tue May 30 10:51:10 2017

@author: noysena
"""
import time
import os
from astropy.io import fits
#from astropy.wcs import WCS
#SEX parameter in dictionary; this will set for diferent TAROT image

# Default configuration file for SExtractor 2.5.0
# EB 2006-07-14
#
# =============================================================================
#-------------------------------- Default general -----------------------------
# =============================================================================
Default_SEX={'CATALOG_NAME' : '/tmp/output.cat', # name of the output catalog
         'CATALOG_TYPE'     : 'ASCII_HEAD',      # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
         'PARAMETERS_NAME'  : '/usr/share/sextractor/tarot_default.param',

#------------------------------- Extraction ----------------------------------
         'DETECT_TYPE'      : 'CCD',    # CCD (linear) or PHOTO (with gamma correction)
         'DETECT_MINAREA'   : '7.0',    # TCA TCH 5.0 TRE 7.0 minimum number of pixels above threshold
         'DETECT_THRESH'    : '5.0',    # TCA TCH 3.0 TRE 5.0 Tuse '1.5' or '23.2,25.0',    # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'ANALYSIS_THRESH'  : '3.0',    # TCA TCH 2.0 TRE 3.0sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'FILTER'           : 'Y',      # apply filter for detection (Y or N)?
         'FILTER_NAME'      : '/usr/share/sextractor/gauss_1.5_3x3.conv',   # name of the file containing the filter
         'DEBLEND_NTHRESH'  : '64',     # Number of deblending sub-thresholds
         'DEBLEND_MINCONT'  : '0.005',  # Minimum contrast parameter for deblending
         'CLEAN'            : 'Y',      # Clean spurious detections? (Y or N)?
         'CLEAN_PARAM'      : '1.0',    # Cleaning efficiency
         'MASK_TYPE'        : 'CORRECT',# type of detection MASKing: can be one of NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
        'PHOT_APERTURES'    : '10',     # MAG_APER aperture diameter(s) in pixels; normal set to 5
        'PHOT_AUTOPARAMS'   : '5., 10.',# MAG_AUTO parameters: <Kron_fact>,<min_radius> '2.5, 3.5'
        'PHOT_PETROPARAMS'  : '5., 10.', # MAG_PETRO parameters: <Petrosian_fact>, <min_radius> '2.5, 3.5'
        'SATUR_LEVEL'       : '55930.0',# level (in ADUs) at which arises saturation
        'MAG_ZEROPOINT'     : '25.0',   # magnitude zero-point
        'MAG_GAMMA'         : '4.0',    # gamma of emulsion (for photographic scans)
        'GAIN'              : '2.0',    # detector gain in e-/ADU
        'PIXEL_SCALE'       : '0',      # Tarot Calern 3.2 size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
        'SEEING_FWHM'       : '1.84',   # stellar FWHM in arcsec
        'STARNNW_NAME'      : '/usr/share/sextractor/default.nnw',    # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------
        'BACK_SIZE'         : '64',     # Background mesh: <size> or <width>,<height>
        'BACK_FILTERSIZE'   : '3',      # Background filter: <size> or <width>,<height>
        'BACKPHOTO_TYPE'    : 'GLOBAL', # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
        'CHECKIMAGE_TYPE'   : 'APERTURES',# can be NONE, BACKGROUND, BACKGROUND_RMS,
                                          # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                          # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                          # or APERTURES
        'CHECKIMAGE_NAME'   : '/tmp/check.fits',     # Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
        'MEMORY_OBJSTACK'   : '6000',     # number of objects in stack
        'MEMORY_PIXSTACK'   : '4194304',  # number of pixels in stack
        'MEMORY_BUFSIZE'    : '2048',     # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
        'VERBOSE_TYPE'      : 'NORMAL',   # can be QUIET, NORMAL or FULL
        'WRITE_XML'         : 'N',        # Write XML file (Y/N)?
        'XML_NAME'          : 'sex.xml'}  # Filename for XML output


# =============================================================================
#-------------------------------- Default 2k-----------------------------------
# =============================================================================
Default_SEX_2k={'CATALOG_NAME' : '/tmp/output.cat', # name of the output catalog
         'CATALOG_TYPE'     : 'ASCII_HEAD',             # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
         'PARAMETERS_NAME'  : '/usr/share/sextractor/tarot_default.param',

#------------------------------- Extraction ----------------------------------
         'DETECT_TYPE'      : 'CCD',  # CCD (linear) or PHOTO (with gamma correction)
         'DETECT_MINAREA'   : '4.0',  # TCA TCH 5.0 TRE 7.0 minimum number of pixels above threshold
         'DETECT_THRESH'    : '2.0',  # TCA TCH 3.0 TRE 5.0 Tuse '1.5' or '23.2,25.0',    # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'ANALYSIS_THRESH'  : '2.0',  # TCA TCH 2.0 TRE 3.0sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'FILTER'           : 'Y',    # apply filter for detection (Y or N)?
         'FILTER_NAME'      : '/usr/share/sextractor/gauss_1.5_3x3.conv',   # name of the file containing the filter
         'DEBLEND_NTHRESH'  : '32',   # Number of deblending sub-thresholds
         'DEBLEND_MINCONT'  : '0.005',# Minimum contrast parameter for deblending
         'CLEAN'            : 'Y',    # Clean spurious detections? (Y or N)?
         'CLEAN_PARAM'      : '1.0',  # Cleaning efficiency
         'MASK_TYPE'        : 'CORRECT', # type of detection MASKing: can be one of NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
        'PHOT_APERTURES'    : '5',       # MAG_APER aperture diameter(s) in pixels; normal set to 5
        'PHOT_AUTOPARAMS'   : '2.5, 3.5',# MAG_AUTO parameters: <Kron_fact>,<min_radius> '2.5, 3.5'
        'PHOT_PETROPARAMS'  : '2.5, 3.5',# MAG_PETRO parameters: <Petrosian_fact>, <min_radius> '2.5, 3.5'
        'SATUR_LEVEL'       : '55930.0', # level (in ADUs) at which arises saturation
        'MAG_ZEROPOINT'     : '25.0',    # magnitude zero-point
        'MAG_GAMMA'         : '4.0',     # gamma of emulsion (for photographic scans)
        'GAIN'              : '2.0',     # detector gain in e-/ADU
        'PIXEL_SCALE'       : '3.35',    # Tarot Calern 3.35, TRE 3.73 size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
        'SEEING_FWHM'       : '1.84',    # stellar FWHM in arcsec
        'STARNNW_NAME'      : '/usr/share/sextractor/default.nnw',    # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------
        'BACK_SIZE'         : '64',      # Background mesh: <size> or <width>,<height>
        'BACK_FILTERSIZE'   : '3',       # Background filter: <size> or <width>,<height>
        'BACKPHOTO_TYPE'    : 'GLOBAL',  # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
        'CHECKIMAGE_TYPE'   : 'APERTURES',# can be NONE, BACKGROUND, BACKGROUND_RMS,
                                          # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                          # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                          # or APERTURES
        'CHECKIMAGE_NAME'   : '/tmp/check.fits',# Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
        'MEMORY_OBJSTACK'   : '6000',     # number of objects in stack
        'MEMORY_PIXSTACK'   : '4194304',  # number of pixels in stack
        'MEMORY_BUFSIZE'    : '2048',     # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
        'VERBOSE_TYPE'      : 'NORMAL',   # can be QUIET, NORMAL or FULL
        'WRITE_XML'         : 'N',        # Write XML file (Y/N)?
        'XML_NAME'          : 'sex.xml'}  # Filename for XML output

# =============================================================================
#-------------------------------- Default 4k-----------------------------------
# =============================================================================

Default_SEX_4k={'CATALOG_NAME' : '/tmp/output.cat', # name of the output catalog
         'CATALOG_TYPE'     : 'ASCII_HEAD',# NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
         'PARAMETERS_NAME'  : '/usr/share/sextractor/tarot_default.param',

#------------------------------- Extraction ----------------------------------
         'DETECT_TYPE'      : 'CCD',  # CCD (linear) or PHOTO (with gamma correction)
         'DETECT_MINAREA'   : '4.0',  # TRE 7.0 minimum number of pixels above threshold
         'DETECT_THRESH'    : '2.0',  # TRE 5.0 Tuse '1.5' or '23.2,25.0',    # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'ANALYSIS_THRESH'  : '2.0',  # TRE 3.0sigmas> or <threshold>,<ZP> in mag.arcsec-2
         'FILTER'           : 'Y',    # apply filter for detection (Y or N)?
         'FILTER_NAME'      : '/usr/share/sextractor/gauss_1.5_3x3.conv',   # name of the file containing the filter
         'DEBLEND_NTHRESH'  : '64',   # Number of deblending sub-thresholds
         'DEBLEND_MINCONT'  : '0.005',# Minimum contrast parameter for deblending
         'CLEAN'            : 'Y',    # Clean spurious detections? (Y or N)?
         'CLEAN_PARAM'      : '1.0',  # Cleaning efficiency
         'MASK_TYPE'        : 'CORRECT', # type of detection MASKing: can be one of NONE, BLANK or CORRECT
 
#------------------------------ Photometry -----------------------------------
        'PHOT_APERTURES'    : '10',      # MAG_APER aperture diameter(s) in pixels; normal set to 5
        'PHOT_AUTOPARAMS'   : '4.5, 5.5',# MAG_AUTO parameters: <Kron_fact>,<min_radius> '2.5, 3.5'
        'PHOT_PETROPARAMS'  : '4.5, 5.5',# MAG_PETRO parameters: <Petrosian_fact>, <min_radius> '2.5, 3.5'
        'SATUR_LEVEL'       : '55930.0', # level (in ADUs) at which arises saturation
        'MAG_ZEROPOINT'     : '25.0',    # magnitude zero-point
        'MAG_GAMMA'         : '4.0',     # gamma of emulsion (for photographic scans)
        'GAIN'              : '2.0',     # detector gain in e-/ADU
        'PIXEL_SCALE'       : '3.73',    # Tarot Calern 3.35, TRE 3.73 size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
        'SEEING_FWHM'       : '2.0',     # stellar FWHM in arcsec
        'STARNNW_NAME'      : '/usr/share/sextractor/default.nnw',    # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------
        'BACK_SIZE'         : '64',      # Background mesh: <size> or <width>,<height>
        'BACK_FILTERSIZE'   : '3',       # Background filter: <size> or <width>,<height>
        'BACKPHOTO_TYPE'    : 'GLOBAL',  # can be GLOBAL or LOCAL
 
#------------------------------ Check Image ----------------------------------
        'CHECKIMAGE_TYPE'   : 'APERTURES',# can be NONE, BACKGROUND, BACKGROUND_RMS,
                                          # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                          # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                          # or APERTURES
        'CHECKIMAGE_NAME'   : '/tmp/check.fits',# Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
        'MEMORY_OBJSTACK'   : '6000',     # number of objects in stack
        'MEMORY_PIXSTACK'   : '4194304',  # number of pixels in stack
        'MEMORY_BUFSIZE'    : '4096',     # number of lines in buffer
 
#----------------------------- Miscellaneous ---------------------------------
        'VERBOSE_TYPE'      : 'NORMAL',   # can be QUIET, NORMAL or FULL
        'WRITE_XML'         : 'N',        # Write XML file (Y/N)?
        'XML_NAME'          : 'sex.xml'}  # Filename for XML output

# =============================================================================
# Function to read and export default sex
# =============================================================================
#print SEXpara.keys(), SEXpara.values()
#built parameter.sex for SExtractor and export to /tmp
def CHECK_STAT(input_value, output_value):
    input_value = output_value
    return input_value

def DEFAULT_SEX(filefits):
    hdu = fits.open(filefits)
    
    try:
        CHECK_STAT(Default_SEX['SEEING_FWHM'], hdu[0].header['FWHM'])      #Star/Galaxy separation
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX['SATUR_LEVEL'], hdu[0].header['DATAMAX'])
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX['CHECKIMAGE_NAME'], hdu[0].header['FILENAME'])
    except(KeyError):
        pass
    
    Sex = open('/tmp/tarot_default.sex', 'w')
    for i, j in zip(Default_SEX.keys(), Default_SEX.values()):
        DSEX = ("%s %s\n" %('{0:<20}'.format(i), '{0:<50}'.format(j)))
        Sex.write(DSEX)
    Sex.close()
    while not os.path.isfile('/tmp/tarot_default.sex'):
        time.sleep(1)

    return Default_SEX

def DEFAULT_SEX_2k(filefits):
    hdu = fits.open(filefits)
    
    try:
        CHECK_STAT(Default_SEX_2k['SEEING_FWHM'], hdu[0].header['FWHM'])      #Star/Galaxy separation
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX_2k['SATUR_LEVEL'], hdu[0].header['DATAMAX'])
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX_2k['CHECKIMAGE_NAME'], hdu[0].header['FILENAME'])
    except(KeyError):
        pass
    
    Sex = open('/tmp/tarot_default.sex', 'w')
    for i, j in zip(Default_SEX_2k.keys(), Default_SEX_2k.values()):
        DSEX = ("%s %s\n" %('{0:<20}'.format(i), '{0:<50}'.format(j)))
        Sex.write(DSEX)
    Sex.close()
    while not os.path.isfile('/tmp/tarot_default.sex'):
        time.sleep(1)

    return Default_SEX

def DEFAULT_SEX_4k(filefits):
    hdu = fits.open(filefits)
    try:
        CHECK_STAT(Default_SEX_4k['SEEING_FWHM'], hdu[0].header['FWHM'])      #Star/Galaxy separation
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX_4k['SATUR_LEVEL'], hdu[0].header['DATAMAX'])
    except(KeyError):
        pass
    
    try:
        CHECK_STAT(Default_SEX_4k['CHECKIMAGE_NAME'], hdu[0].header['FILENAME'])
    except(KeyError):
        pass
    
    Sex = open('/tmp/tarot_default.sex', 'w')
    for i, j in zip(Default_SEX_4k.keys(), Default_SEX_4k.values()):
        DSEX = ("%s %s\n" %('{0:<20}'.format(i), '{0:<50}'.format(j)))
        Sex.write(DSEX)
    Sex.close()
    while not os.path.isfile('/tmp/tarot_default.sex'):
        time.sleep(1)

    return Default_SEX
