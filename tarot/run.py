#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:10:24 2018

@author: kanthanakorn
"""
# =============================================================================
# How to run
# import tarot
# print(tarot.TAROT_PIP.MAIN())
# =============================================================================

# =============================================================================
# Import from Astropy packages and Python packages
# =============================================================================
import os
import sys
import time
import copy
import timeit
import subprocess
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.table import Table
from astropy.table import Column
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky

# =============================================================================
# Import from TSP
# =============================================================================
from .utility.test_run import Test_Run
from .filter.iteration import Iterate
from .catalog.LGaia2_Tarot10 import LoadCata
from .photometry.star_ref import Mag_Ref_Table
from .filter.duplication import Duplicate_Match
from .catalog.LGlade_Tarot10 import LoadCataGalaxy
from .utility.limit_in_fits import Limit_Coor_to_Fits
from .filter.check_unknown_source import Check_Object
from .filter.check_unknown_source import Check_Object_NOMAD

speed = False
if speed:
    from .sextractor.Default_SEx_for_TAROT_Fast_Search import DEFAULT_SEX
    from .sextractor.Default_SEx_for_TAROT_Fast_Search import DEFAULT_SEX_2k
    from .sextractor.Default_SEx_for_TAROT_Fast_Search import DEFAULT_SEX_4k
else:
    from .sextractor.Default_SEx_for_TAROT import DEFAULT_SEX
    from .sextractor.Default_SEx_for_TAROT import DEFAULT_SEX_2k
    from .sextractor.Default_SEx_for_TAROT import DEFAULT_SEX_4k

class TAROT_PIP(object):
#    Class and Main() is help to call packpage easier with . (search in this and after folder) and .. (go up to upper folder and search for package)
    def MAIN():
        return("Transient search package for searching an unkown candidate in TAROT images")
        
        """ TAROT pipeline cadidate searching attributes;
        fitsfile : import FITS
        dir_fitsfile    : directory of FITS
        name            : name of FITS
        fname           : prefix name with out type
        catalog_name    : referent catalog
        SExoutput_name  : Source Extractor file
        candidate_name  : candidate file
        SExdata         : Source Extractor data
        catalog         : catalog data
        tbmagnitude     : table of magnitude
        coordinate_data : coordinate data
        coordinate_catalog : coordinate catalog data
        HDU_header      : header FITS
        wcs             : WCS FITS
        fitsnew         : FITS after distortion correction
        new_HDU_header  : header FITS after distortion correction
        new_wcs         : WCS FITS after distortion correction
        idx             : candidate index after cross-matching
        d2d             : angular distance obtaining from cross-matching
        d3d             : souce distance if we have distance of our data
        matches         : catalog that matches with data
        idx_candi_data  : idx of unknow candidate
        image           : image quality confirmation
        """
# =============================================================================
# Funtion running in TSP
# =============================================================================
    def __init__(self, fitsfile):
        self.fitsfile = fitsfile
        self.dir_fitsfile = os.path.dirname(self.fitsfile)
        self.name = os.path.basename(self.fitsfile)
        self.fname, self.lname = self.name.rsplit(".")
        self.catalog_name = os.path.join(self.dir_fitsfile, self.fname +'.GaiaDR2.dat')
        self.SExoutput_name = os.path.join(self.dir_fitsfile, self.fname +'.cat')
        self.candidate_name = os.path.join(self.dir_fitsfile, self.fname +'.candi.dat')
        self.candidate_html = os.path.join(self.dir_fitsfile, self.fname +'.html')
        self.fitsnew = os.path.join(self.dir_fitsfile, self.fname+'.new')
        self.test = False
        
        self.HDU_header = ''
        self.new_HDU_header = ''
        self.wcs = ''
        self.new_wcs = ''
        
        self.catalog = ''
        self.SExdata = ''
        self.coordinate_catalog = ''
        self.catalog_galaxy_table = ''
        
        
        self.idx = ''
        self.d2d = ''
        self.d3d = ''
        self.matches = ''
        self.Cdata = ''
        self.idx_candi_data = ''
        
        self.confirmed_candidate = ''
        self.confirmed_candidate_idx = ''
        self.d2d_candidate = ''
        
        
    def Check_Test(self):
        self.test = Test_Run()
    
        print("Test status %s\nRun on machine: %s" %(self.test, os.getcwd()))
        print(self.fitsfile)

    """Read FITS"""
    def Fits_Info(self):
        hdu = fits.open(self.fitsfile)
        self.HDU_header = hdu[0].header
        self.wcs = WCS(hdu[0].header)
        hdu.close()
    
    """Read FITS after distortion correction with AStrometry in extention NEW"""
    def Fits_New_Info(self):
        hdu = fits.open(self.fitsnew)
        self.new_HDU_header = hdu[0].header
        self.new_wcs = WCS(hdu[0].header)
        hdu.close()

    """Run Astrometry at tweak 5th to get FITS in extention NEW"""
    def Run_Astrometry(self):
        """Astrometry is run by bash
        clean *.axy and *.wcs, !!!--> lead to solved error, keep axy, wcs"""
#        Xpixel = self.HDU_header['NAXIS1']
#        Ypixel = self.HDU_header['NAXIS2']
        
#        RA0, DEC0 = self.wcs.all_pix2world(0, 0, 0)
#        RA, DEC = self.wcs.all_pix2world(Xpixel/2.0, Ypixel/2.0, 0)

        """Calculate radius of FOV for catalog frame"""
#        Corner_frame = SkyCoord(RA0, DEC0, frame = 'fk5', unit = (u.deg, u.deg))
#        Center_frame = SkyCoord(RA, DEC, frame = 'fk5', unit = (u.deg, u.deg))
#        FOV = Corner_frame.separation(Center_frame)
#        Radius = FOV.deg + 0.10
        Radius = 6.0
        
        a=0
        try:
            TAROT_PIP.Fits_New_Info(self)
        except:
#            runcode = ("solve-field --fits-image --guess-scale --no-plots --skip-solved --corr none --index-xyls none --match none --solved none --rdls none --overwrite --tweak-order 5 --ra %s --dec %s --radius %s --dir %s %s" 
            runcode = ("solve-field --fits-image --guess-scale --no-plots\
                       --skip-solved --corr none --index-xyls none --match none\
                       --solved none --rdls none --overwrite --tweak-order 5 \
                       --ra %s --dec %s --radius %s --dir %s %s" 
                       %(self.HDU_header['RA'], self.HDU_header['DEC'],
                       Radius, self.dir_fitsfile, self.fitsfile))
            subprocess.Popen(runcode, stdout=subprocess.PIPE, shell=True)
        
        while not os.path.exists(self.fitsnew):
            a += 1
            print("Field is solving...%s" %(a), end='\r')
            time.sleep(1)
            
            if a == 360:
                with open(self.candidate_name, 'w') as bad_cadi:
                    bad_cadi.write("pbcalibration")
#                    os.remove(self.fitsfile) #we have loop to unzip, hash out
                break
                
            if os.path.exists(self.fitsnew):
                TAROT_PIP.Fits_New_Info(self)
                break
            else:
                pass
        else:
            pass
            
    """Run SExtractor to get XY pixel and flux of sources in TAROT image"""
    def Run_SExtractor(self):
        
        Xpixel = self.HDU_header['NAXIS1']
        Ypixel = self.HDU_header['NAXIS2']
        
        Xlim = Xpixel - 60
        Ylim = Ypixel - 60
        
        print("Check 2048", Xpixel == 2048)
        print("Check 4096", Xpixel == 4096)
        
        if Xpixel == 4096:
            DEFAULT_SEX_4k(self.fitsfile)
        elif Xpixel == 2048:
            DEFAULT_SEX_2k(self.fitsfile)
        else:
            DEFAULT_SEX(self.fitsfile)
            
        
        if os.path.isfile('/tmp/output.cat'):
            os.remove('/tmp/output.cat')
        else:
            pass
            
        runcode = ("sextractor -c /tmp/tarot_default.sex " + self.fitsfile)
        run = subprocess.Popen(runcode, stdout=subprocess.PIPE, shell=True)
        out, err = run.communicate()
        
        """Use while not to check that is file existant """
        if not os.path.exists("/tmp/output.cat"):
            time.sleep(5)
        else:
            pass
            
        move_cat = ('cp /tmp/output.cat %s' %self.SExoutput_name)
        run = subprocess.Popen(move_cat, stdout=subprocess.PIPE, shell = True)
        out, err = run.communicate()

        try:
            self.SExdata = Table.read(self.SExoutput_name,
                                      format='ascii.sextractor')
            
            xlimit_left = np.where(self.SExdata["X_IMAGE"] < 60)
            self.SExdata.remove_rows(xlimit_left)
            xlimit_right = np.where(self.SExdata["X_IMAGE"] > Xlim)
            self.SExdata.remove_rows(xlimit_right)
            
            ylimit_down = np.where(self.SExdata["Y_IMAGE"] < 60)
            self.SExdata.remove_rows(ylimit_down)
            ylimit_up = np.where(self.SExdata["Y_IMAGE"] > Ylim)
            self.SExdata.remove_rows(ylimit_up)
            
        except:
            print("SExtractor error")

    def Load_Galaxies(self, RA, DEC, Radius):
        galaxy_table = LoadCataGalaxy(RA, DEC, round(Radius, 2))
        filter_null = galaxy_table['dist'] != 'null'
        self.catalog_galaxy_table = galaxy_table[filter_null]
        del galaxy_table
        
        keep_gal = []
        if True:
            for i, j in enumerate(self.catalog_galaxy_table['dist']):
                if float(j) < 150:
                    keep_gal.append(i)
            self.catalog_galaxy_table = self.catalog_galaxy_table[keep_gal]
            
#        print(self.catalog_galaxy_table)
        
        self.catalog_galaxy_table = Limit_Coor_to_Fits(filefits = self.fitsfile,
                                                       catalog = self.catalog_galaxy_table,
                                                       ratext = 'RA',
                                                       detext = 'dec')
        return self.catalog_galaxy_table

    def Load_Catalog(self):
        TAROT_PIP.Check_Test(self)
        #No new file, exit
        if not os.path.exists(self.fitsnew):
            sys.exit()
        print("Recieving catalog from Gaia DR2 ...")
        
        Xpixel = self.new_HDU_header['NAXIS1']
        Ypixel = self.new_HDU_header['NAXIS2']
        
        RA0, DEC0 = self.new_wcs.all_pix2world(0, 0, 0)
        RA, DEC = self.new_wcs.all_pix2world(Xpixel/2.0, Ypixel/2.0, 0)

        """Calculate radius of FOV for catalog frame"""
        Corner_frame = SkyCoord(RA0, DEC0, frame = 'fk5', unit = (u.deg, u.deg))
        Center_frame = SkyCoord(RA, DEC, frame = 'fk5', unit = (u.deg, u.deg))
        FOV = Corner_frame.separation(Center_frame)
        Radius = FOV.deg + 0.10
        print ("FOV of catalog RA DEC Radius: %s %s %s"
               %(np.round(RA, 2),
                 np.round(DEC, 2),
                 round(Radius, 2)))
        
        if self.test:
# =============================================================================
#             Test catalog in Noysena machine
# =============================================================================
            myhome = os.path.expanduser('~')
            fpath = os.path.join(myhome, "/home/kanthanakorn/Documents/Data/Test/")
            #TRE
#            f = os.path.join(fpath,
#                             "IM_20170105_165445540_000000_44954708_cata.csv")
            #TCH
            f = os.path.join(fpath,
                             "IM_20191207_023356000_000002_26943902.cata.csv")
            
            catalog_table=Table.read(f, format = 'csv')
            
            #load galaxies
            self.catalog_galaxy_table = TAROT_PIP.Load_Galaxies(self, RA, DEC, Radius)
            
# =============================================================================
#             Catalog load from link in case we have a very fast download
# =============================================================================
            #Load catalog form xml is slow
#            url_gaia_dr1 = "http://0.0.0.0:1112/cgi-bin/cs.py?RA={:>3.6f}&DEC={:>3.6f}&SR={:>2.3f}".format(RA, DEC, Radius)
#            Retrive_URL = parse_single_table(url_gaia_dr1)
#            catalog_votable = Retrive_URL.array
#            catalog_table = Table(catalog_votable.data)
        else:
            #Load catalog directly form file is fast
            #Find the limite magnitude for catalog frame
            catalog_table = LoadCata(RA, DEC, round(Radius, 2))
            
            #Load Galaxy from Glade catalog
            try:
                self.catalog_galaxy_table = TAROT_PIP.Load_Galaxies(self, RA, DEC, Radius)
            except:
                with open(self.candidate_name, 'w') as bad_cadi:
                    bad_cadi.write("load catalog problem")
        
        """Need to obtimize the limitting magnitude by comparing Rmag"""
        
        if Xpixel == 4096:
            limit_mag = np.where(catalog_table['phot_g_mean_mag'] <= 17.0)
        elif Xpixel == 2048: 
            limit_mag = np.where(catalog_table['phot_g_mean_mag'] <= 18.0)
        else:
            limit_mag = np.where(catalog_table['phot_g_mean_mag'] <= 19.0)
            
        self.catalog = copy.copy(catalog_table[limit_mag])
        catalog_table = ""
        print("{:>d} rows of Catalog.".format(len(self.catalog)))

    """Read data and catalog andmake skycoordinate"""
    def Read_Data(self):
#        self.SExdata = Table.read(self.SExoutput_name,
#                                  format='ascii.sextractor')
        
        TAROT_PIP.Load_Catalog(self)
        
        self.coordinate_catalog = SkyCoord(self.catalog['ra'],
                                           self.catalog['dec'],
                                           frame = 'icrs',
                                           unit = (u.deg, u.deg))
        
        #Load data and built sky coordinate with WCS of 5th SIP
        RA, DEC = self.new_wcs.all_pix2world(self.SExdata['X_IMAGE'],
                                             self.SExdata['Y_IMAGE'], 0)
        fk5_coordinate_data = SkyCoord(np.round(RA, 3),
                                       np.round(DEC,3),
                                       frame = 'fk5',
                                       unit = (u.deg, u.deg))
        self.coordinate_data = fk5_coordinate_data.transform_to('icrs')
        
        self.SExdata["X_WORLD"] = np.round(self.coordinate_data.ra.deg, 3)
        self.SExdata["Y_WORLD"] = np.round(self.coordinate_data.dec.deg, 3)
        
        #Convert to 'MAG_AUTO' to relative mag using star near 'RA', 'DEC' in hdu
        try:
            self.SExdata = Mag_Ref_Table(self.fitsnew, self.SExdata)
        except:
            print("Error convert magnitude of Mag_Ref_Table")

    # I don't need it before iteration as I Matching in Iterate script
    def Data_Match_Catalog(self, show_cpu_time=False):
        if show_cpu_time:
            start = timeit.default_timer()
            TAROT_PIP.Read_Data(self)
            self.idx,self.d2d,self.d3d = match_coordinates_sky(self.coordinate_data,
                                                               self.coordinate_catalog)
            self.matches = self.coordinate_catalog[self.idx]
            stop = timeit.default_timer();
            print("Matching had done in {:> 2.2f} seconds".format(stop-start))
            print("Angular separation : {:> 2.2f} arcsec,\
                  STD: {:> 2.2f}".format(np.median(self.d2d.arcsec),
                  np.std(self.d2d.arcsec)))
        else:
            TAROT_PIP.Read_Data(self)
            self.idx,self.d2d,self.d3d = match_coordinates_sky(self.coordinate_data,
                                                               self.coordinate_catalog)
            self.matches = self.coordinate_catalog[self.idx]

    def Iteration(self, show_cpu_time=False):
        if show_cpu_time:
            start = timeit.default_timer()
            TAROT_PIP.Read_Data(self)
            self.idx, self.d2d, self.d3d,self.matches, self.Cdata = Iterate(self.coordinate_data,
                                                                            self.coordinate_catalog)
            stop = timeit.default_timer();
            print("Matching had done in {:> 2.2f} seconds".format(stop-start))
            print("Angular separation: {:> 2.2f} arcsec, STD: {:> 2.2f}".format(np.median(self.d2d.arcsec),
                  np.std(self.d2d.arcsec)))
        else:
            TAROT_PIP.Read_Data(self)
            self.idx, self.d2d, self.d3d,self.matches, self.Cdata = Iterate(self.coordinate_data,
                                                                        self.coordinate_catalog)
        
            stop = timeit.default_timer();
        
    def Search_New_Object(self):
        """
        Input : Standard derivation number, indixes, angular distance, data
        Output: Indices of catalog, indixes of data
        """
        sigma = 3
        missmatch = np.median(self.d2d.arcsec) + sigma*self.d2d.arcsec.std()
        print("Missmatch value: %s" %missmatch)
        idx_candi_data_0 = np.where(self.d2d.arcsec >= missmatch)
        #data indices of data due to it is created from np.where
        self.idx_candi_data = np.array(idx_candi_data_0).flatten()

        if False:
            print("No. of unknow source: %s" %str(len(self.idx_candi_data)))
        else:
            pass
        
        print("Unknown source over median: %s" %str(len(self.idx_candi_data)))

    def Add_Duplicate_Object(self):
        duplicate_idx = Duplicate_Match(self.idx)
        print("Number of duplicate source %s" %len(duplicate_idx))
        # Change to True to see number of duplication
        if False:
            print(duplicate_idx)
        else:
            pass
        
        if len(duplicate_idx) > 0:
            for idx in duplicate_idx:
                if idx in self.idx_candi_data:
                    pass
                else:
                    self.idx_candi_data = np.append(self.idx_candi_data,(idx))
        else:
            pass
        #If need List, use sorted(), result in []
        self.idx_candi_data = np.sort(self.idx_candi_data)
        print("Number of candidate: %s" %len(self.idx_candi_data))
    
    def Check_Candidate(self):
        Search_by_Nomad = True
        if Search_by_Nomad:
            self.confirmed_candidate_idx = Check_Object_NOMAD(self.SExdata[self.idx_candi_data],
                                                              self.Cdata[self.idx_candi_data])
        else:
            self.confirmed_candidate_idx = Check_Object(self.SExdata[self.idx_candi_data],
                                                              
                                                    self.Cdata[self.idx_candi_data])
        
        self.confirmed_candidate = self.SExdata[self.idx_candi_data[self.confirmed_candidate_idx]]
        self.d2d_candidate = self.d2d.arcsec[self.idx_candi_data[self.confirmed_candidate_idx]]
        
        d2d_candidate_column = Column(np.round(self.d2d_candidate), name = 'AngSep')
        self.confirmed_candidate.add_column(d2d_candidate_column,
                                            index = None)
#        self.confirmed_candidate.write("Candidate.xml",
#                                       format="votable", 
#                                       overwrite=True)
        
        #Eliminate bright object from candidate
        bright_obs = np.where(self.confirmed_candidate['MAG_AUTO'] < 12.)
        self.confirmed_candidate.remove_rows(bright_obs)
        
        self.confirmed_candidate.write(self.candidate_name, 
                                       format="ascii", 
                                       overwrite=True)
        self.confirmed_candidate.write(self.candidate_html, 
                                       format="ascii", 
                                       overwrite=True)
        
        print(self.confirmed_candidate)
        print("No. of possible candidate: ", len(self.confirmed_candidate))
        print("ID: ", self.confirmed_candidate["NUMBER"])
        
