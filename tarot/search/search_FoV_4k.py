#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 12:29:54 2019

@author: kanthanakorn
"""
import copy
import time
import os
import sys
import glob
import numpy as np
from numba import jit

from skimage import exposure
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.table import Table, Column, vstack
from tarot.utility.sub_fits import SUB_FITS
from tarot.filter.ml_source import ML_SOURCE
from tarot.database.TSP_SQL import Insert_Image
from tarot.utility.general_funtion import sub_image
from tarot.utility.general_funtion import complet_array
from tarot.catalog.load_image_from_DSS import LOAD_DSS_IMAGE
from tarot.utility.general_funtion import Plot_Image_Candi_Solar_Obj


info = """ Beware, if there is not collect FoV or enough FoV as mention
           in Fits name, then error occur """
           
class TAROT_FoV(object):
    @jit(nopython=True)
    def MAIN():
        return("Searching FoV in TAROT images")
    # Define paramater and value
    def __init__(self, Path):
        self.path = Path
        self.fits = None
        self.fitslist = None
        self.HDU_header = None
        self.sname = ''
        self.wcs = None
        self.FoV = None
        self.table = Table()
        self.max_FoV = 0
        self.list_FoV = [] # All FoV in each image
        self.FoVs = None   # Take form kmeans after group from many FoVs
        self.path_PNG_list = "/tmp/keep_PNG"
        self.path_PNG_folder= "/tmp/keep_PNG"
#        self.clf =  ML_SOURCE('source_identification.joblib')
        
    # Search all FITS in folder
    def List_File(self):
        if os.access(self.path, os.R_OK): #Ture
            findit = os.path.join(self.path, '*.new')
            subfit = os.path.join(self.path, '*_sub.fits')
            # fitslist = []
            fitslist = sorted(glob.glob(findit))
            subflist = sorted(glob.glob(subfit))
            if len(fitslist) > 0:
            #Remove subfits that obtain from DSS/subimage
                for sf in subflist:
                    if sf in fitslist:
                        fitslist.remove(sf)
                    else:
                        pass
       
        if len(fitslist) > 0:
            for fit in fitslist:
                hdu = fits.open(fit)
                hdu_header = hdu[0].header
                hdu.close()
                
                try:
                    hdu_header["USER"]
                    if hdu_header["USER"] == 'alert_grandma':
                        pass
                    else:
                        fitslist.remove(fit)
                except(KeyError):
                    fitslist.remove(fit)
            #Remove subfits that obtain from DSS/subimage

        else:
            fitslist = None
        self.fitslist = fitslist
        
    # Read FITS
    def Fits_Info(self):
        with fits.open(self.fits) as hdu:
            self.HDU_header = hdu[0].header
            self.wcs = WCS(hdu[0].header)
            self.data = hdu[0].data
            self.sname = self.HDU_header['SNAME']
            
            ra = round(self.HDU_header['CRVAL1'],4)
            dec = round(self.HDU_header['CRVAL2'],4)
            self.FoV = [ra, dec]
    
    #Read FoV form SNAME in hearder
    def Read_FoV(self):
        try:
            #update checking tiles 20July2019
            tile = self.sname.split("_")[-1]
            # TRE header update 25/04/2019
#            tile = self.sname.split("_")[3]
            
            # TRE old header
#            fov = self.sname.split("_")[1]
        except(ValueError, IndexError):
            print("Read FoV error: %s " %self.fits)
            sys.exit()
        return int(tile)
  
    # Output FoV table
    def FoV_Table(self):
        hdu_name = ["sname", "FILENAME", "Time-OBS", "RA", "DEC", "COORDSYS",
                    "FWHM", "EXPOSURE", "BGMEAN", "BGSIGMA", "FILTER", "TELESCOP"]

    # Read hud and check Number of FoV then list out FoV table
        Max_FoV_TCA = 0
        Max_FoV_TCH = 0
        Max_FoV_TRE = 0
        
        # Call fitslist of fits in folder
        TAROT_FoV.List_File(self)
        
        for i, j in enumerate(self.fitslist):
            self.fits = j
#            self.HDU_header, self.wcs, self.FoV = 
            TAROT_FoV.Fits_Info(self)
            self.list_FoV.append(self.FoV)
    
    # Check hdu TRE
            try:
                tele = self.HDU_header["TELESCOP"]
            except(KeyError):
                tele = "TAROT"
        
            if tele == "TAROT CALERN":
                FoV_TCA = TAROT_FoV.Read_FoV(self)
        
        #Check FoV, it is maximum for TCA or not
                if FoV_TCA > Max_FoV_TCA:
                    Max_FoV_TCA = FoV_TCA
                else:
                    pass
        
        # Check hdu TCH
            elif tele == "TAROT CHILI":
                FoV_TCH = TAROT_FoV.Read_FoV(self)
        
        # Check FoV, it is maximum for TCH or not
                if FoV_TCH > Max_FoV_TCH:
                    Max_FoV_TCH = FoV_TCH
                else:
                    pass
            
        # Check hdu TRE
            elif tele == "TAROT REUNION":
                FoV_TRE = TAROT_FoV.Read_FoV(self)

        # Check Fov, it is maximum for TRE or not
                if FoV_TRE > Max_FoV_TRE:
                    Max_FoV_TRE = FoV_TRE
                else:
                    pass
            else:
                pass
      # Read hud info in each fits for FoV table      
            try:
                fwhm = self.HDU_header["FWHM"]
            except(KeyError):
                fwhm = 99 # NO info for fwhm then 99
                
     # Read Coordinate system --> we chage FKR to ICRS in TAROT.run.Read_Data
     # Here we read it directly, so it must be the same FoV
#            try:
#                self.HDU_header["COORDSYS"] = "ICRS"
#            except(KeyError):
#                pass
            
            try:
                fillter = self.HDU_header["FILTER"]
                if fillter == "N":
                    self.HDU_header["FILTER"] = "Clear"
                elif fillter == "C":
                    self.HDU_header["FILTER"] = "Clear"
                elif fillter == "NoFilter":
                    self.HDU_header["FILTER"] = "Clear"
                else:
                    pass
            except(KeyError):
                self.HDU_header["FILTER"] = "Clear"
                
            try:
                self.HDU_header["COORDSYS"]
            except(KeyError):
                self.HDU_header["COORDSYS"] = 'fk5'
                
       # Build      
            hdu_info = np.array([self.HDU_header["SNAME"],
                                 self.HDU_header["FILENAME"],
                                 self.HDU_header["DATE-OBS"][11:23],
                                 np.around(self.HDU_header["RA"], decimals = 3),
                                 np.around(self.HDU_header["DEC"], decimals = 3),
                                 self.HDU_header["COORDSYS"],
                                 fwhm,
                                 self.HDU_header["EXPOSURE"],
                                 '{:0.3f}'.format(self.HDU_header["BGMEAN"]),
                                 '{:0.3f}'.format(self.HDU_header["BGSIGMA"]),
                                 self.HDU_header["FILTER"], tele])

            if i == 0:
                self.table = Table(hdu_info, names=hdu_name,
                                   masked=True,
                                   meta={'name':'GW counterpart image'},
                                   dtype=('S18', 'S48', 'S12', 'f8', 'f8', 'S4',
                                          'f8', 'f8', 'f8', 'f8', 'S6', 'S13'))
            else:
                self.table.add_row(hdu_info)
        
        # Sum number of FoV from TCA, TCH, TRE to get maximum 
        self.max_FoV = Max_FoV_TCA + Max_FoV_TCH + Max_FoV_TRE
   
        self.table["RA"].unit = u.deg
        self.table["DEC"].unit = u.deg
        self.table["FWHM"].unit = u.arcsec
        self.table["EXPOSURE"].unit = u.s
        
#        return self.table, self.max_FoV, np.array(self.list_FoV)
    
    # Predict and group FoV in folder
    def Predict_FoV(self):
#        self.table, self.max_FoV, FoV_data = 
        TAROT_FoV.FoV_Table(self)
# =============================================================================
#         Predict FoV with SNAME
# =============================================================================
        keep_sn = []
        for sn in self.table:
            """ Change also in line 205 """
            # TRE header update 20/07/2019
            keep_sn.append(int(sn["sname"].split("_")[-1]))
            # TRE header update 25/04/2019
#            keep_sn.append(int(sn["sname"].split("_")[3]))
            
            # TRE old header
#            keep_sn.append(int(sn["sname"].split("_")[1]))
        
        self.FoVs = keep_sn
        col_search_FoV = Column(keep_sn, name="FoV", dtype=int)
        self.table.add_column(col_search_FoV, index=0)

# =============================================================================
#     It call TAROT_FoV.Predict_FoV(self) and call again when recall fn
# =============================================================================
    def Add_Candidate(self):
        TAROT_FoV.Predict_FoV(self)
        app_candi_list = []
        for i in self.table["FILENAME"]:
            fname = i.replace('.fit', '.candi.dat')
            candi_list = os.path.join(self.path, fname)
            if os.path.isfile(candi_list):
                read_candi_list = Table.read(candi_list, format='ascii')
                app_candi_list.append(len(read_candi_list))
            else:
                app_candi_list.append(0) # Either failed search or non candi
                    
        col_candi = Column(app_candi_list, name="Candidate", dtype=int)
        self.table.add_column(col_candi) #add last column is default.
        
    def Export_FoV_Candi(self):
        TAROT_FoV.Add_Candidate(self)
        a = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                      1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        dummy_table = Table(a, names=('NUMBER','FLUX_AUTO','FLUXERR_AUTO',
                            'MAG_AUTO','MAGERR_AUTO','FLUX_BEST','FLUXERR_BEST',
                            'MAG_BEST','MAGERR_BEST','KRON_RADIUS','BACKGROUND',
                            'X_IMAGE','Y_IMAGE','X_WORLD', 'Y_WORLD', 'CXY_IMAGE',
                            'A_IMAGE','B_IMAGE','A_WORLD','B_WORLD','THETA_IMAGE',
                            'ELLIPTICITY','ERRXY_IMAGE','FWHM_IMAGE','FLAGS',
                            'CLASS_STAR','phot_g_mean_mag','ra','dec','Separation'),
                    meta={'name': 'Table Candidate'}, dtype=['S', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8',
                         'f8', 'f8', 'f8', 'f8', 'f8',])
        
# =============================================================================
#         Folder to put PNG and other data in self. must edit for other machine 
# =============================================================================
        if not os.path.exists(self.path_PNG_folder):
            os.makedirs(self.path_PNG_folder)
        else:
            pass
        
        folder_PNG = os.path.join(self.path_PNG_folder, "folder_PNG.txt")            
        if os.path.exists(folder_PNG):
            os.remove(folder_PNG)
        else:
            pass
        
# =============================================================================
#         run from FoV 1 to max + 1 (will exceed to our FoV)
# =============================================================================
        for i in range(1, (max(self.FoVs) + 1)):
            # Create folder to keep candidate
            newpath = os.path.join(self.path, str(i))
            
# =============================================================================
#             Write path for folder contained PNG for using in website.
#             Then I copy it to folder website in website/images/png
# =============================================================================
            with open(folder_PNG, 'a') as folder_i:
                folder_i.write(newpath)
                folder_i.write("\n")
            
            """ define name for ascii data and html for website """
            candi_fov = os.path.join(newpath,
                                     ("all_candidate_FoV_%s.candi.dat" %str(i)))
            candi_fov_html = os.path.join(newpath,
                                          ("all_candidate_FoV_%s.html" %str(i)))
            
            if not os.path.exists(newpath):
                os.makedirs(newpath)
            else:
                print("Path exist: %s\n" %newpath)
                time.sleep(0.5)
            
            same_fov_index = np.where(self.table["FoV"] == i)
            file_list = self.table[same_fov_index]
            
            """ Define name of ascii data and html for website """
            file_list_fov = os.path.join(newpath,
                                         ("table_same_FoV_%s.candi.dat" %str(i)))
            
            file_list_fov_html = os.path.join(newpath,
                                              ("table_same_FoV_%s.html" %str(i)))
            
# =============================================================================
#        Here we include all candiate but we just check duplicate in photometry
# =============================================================================
            for j in range(len(file_list)):
                file = file_list[j]["FILENAME"]
                fname = file.replace('.fit', '.candi.dat')
                file_dat = os.path.join(self.path, fname)
                
# =============================================================================
#                 Read candi.dat file, and build table for vstack
# =============================================================================
                if j == 0:
                    if os.path.exists(file_dat):
                        first_table = Table.read(file_dat,
                                                 format = 'ascii')
                        
                    else:
                        first_table = copy.copy(dummy_table)
                        first_table.remove_row(0)
                else:
                    try:
                        next_table = Table.read(file_dat, format = 'ascii')
                        first_table = vstack([first_table, next_table])
                    except:
                        continue
                
                """ Write table to file """
                first_table.write(candi_fov, format='ascii', overwrite=True)
                first_table.write(candi_fov_html, format='html', overwrite=True)
                
            """ Write file for other script and html for website """
            file_list.write(file_list_fov, format='ascii', overwrite=True)
            file_list.write(file_list_fov_html, format='html', overwrite=True)    
            
# =============================================================================
# ML source    
# =============================================================================
    def ML_Candi(self):
        ml_path = os.path.expanduser('~')
        
        if self.HDU_header == 'TAROT REUNION':
            joblib = os.path.join(ml_path, 'TSP/tarot/filter/source_identification_4k.joblib')
        else:
            joblib = joblib = os.path.join(ml_path, 'TSP/tarot/filter/source_identification_2k.joblib')

        clf = ML_SOURCE(joblib)
            
        for i in range(1, max(self.FoVs)+1):
# =============================================================================
#             Create folder to keep candidate list with PNG info
# =============================================================================
            fov_path = os.path.join(self.path, str(i))
            
            if not os.path.exists(fov_path):
                os.makedirs(fov_path)
            else:
                pass
#                print("Path exist: ", fov_path)
            
            same_fov_index = np.where(self.table["FoV"] == i)
            file_list = self.table[same_fov_index]
            
            # Now exprot PNG
            for j in range(len(file_list)):
                file = file_list[j]["FILENAME"]
                fname = os.path.splitext(file)[0]
                file_dat = os.path.join(self.path, fname) + '.candi.dat'
                file_new = os.path.join(self.path, fname) + '.new'
                
# =============================================================================
#                 Read candi.dat file, and build table for vstack
# =============================================================================
                try:
                    candilist = Table.read(file_dat,
                                             format = 'ascii')
                except:
                    candilist = False
                    pass
                keep = []
                if candilist:
                    for j, k in enumerate(candilist):
                        if os.path.exists(file_new):
                            with fits.open(file_new) as new_hdu:
                                new_data = new_hdu[0].data
#                            print(k["X_IMAGE"], k["Y_IMAGE"])
                            source = sub_image(k["Y_IMAGE"].item(), k["X_IMAGE"].item(), 10, new_data)
                            
                            if clf.predicted == 3:
                                keep.append(j)
                            else:
                                pass
                            if clf.predicted == 0:
                                keep.append(j)
                            else:
                                pass
                            
                            if False:
                                
                                p2, p98 = np.percentile(source, (25, 98))
                                img_contrast = exposure.rescale_intensity(source, in_range=(p2, p98))
                                img = complet_array(img_contrast, 20)
                                clf.What_Source(img)
#                                print(clf.predicted)
                                
                                plt.figure()
                                plt.imshow(img, cmap = 'gray')
                                plt.title("Predicted: %s" %clf.predicted[0])
                                plt.show()
                            else:
                                pass
                        else:
                            pass
                    
                        if len(keep) != 0:
                            candilist.remove_rows(keep)
                            candilist.write(file_dat, format='ascii', overwrite=True)
                        else:
                            pass
                    
                else:
                    pass

# =============================================================================
# Export canidate's image in png        
# =============================================================================
    def Export_PNG_Candi(self):
        path_PNG_list = os.path.join(self.path_PNG_list, "index_PNG_list.txt")
        if os.path.exists(path_PNG_list):
            os.remove(path_PNG_list)
        else:
            pass
        
        for i in range(1, max(self.FoVs)+1):
            # Create folder to keep candidate list with PNG info
            fov_path = os.path.join(self.path, str(i))
            
            if not os.path.exists(fov_path):
                os.makedirs(fov_path)
            else:
                pass
#                print("Path exist: ", fov_path)
            
            same_fov_index = np.where(self.table["FoV"] == i)
            file_list = self.table[same_fov_index]
            
            # Now exprot PNG
            for j in range(len(file_list)):
                file = file_list[j]["FILENAME"]
                fname = os.path.splitext(file)[0]
                file_dat = os.path.join(self.path, fname) + '.candi.dat'
                file_new = os.path.join(self.path, fname) + '.new'
                Insert_Image
# =============================================================================
#                 Read candi.dat file, and build table for vstack
# =============================================================================
                try:
                    candilist = Table.read(file_dat,
                                             format = 'ascii')
                except:
                    candilist = False
                    pass
                
                if candilist:
                    for k in candilist:
                        if os.path.exists(file_new):
# =============================================================================
#                             save PNG in folder and supfits 
# =============================================================================
                            name_PNG = ("%s_%s_%s.PNG" %(fname, str(i), k["NUMBER"].item()))
                            save_PNG_pathfile = os.path.join(fov_path, name_PNG)
                            
                            name_subfit = ("%s_%s_%s.fits" %(fname, str(i), k["NUMBER"].item()))
                            save_subfit = os.path.join(fov_path, name_subfit)
# =============================================================================
#                             Alway write PNG list and distribute for website
# =============================================================================
                            
                            name_PNG_web = ("images/%s/%s_%s_%s.PNG" %(str(i), fname, str(i), k["NUMBER"].item()))
                            name_subfit_web = ("images/%s/%s_%s_%s.fits" %(str(i), fname, str(i), k["NUMBER"].item()))
                            
# =============================================================================
#                             Upload 'png' to database 
# =============================================================================
                            print(self.HDU_header['TELESCOP'])
                            if True:
                                if self.HDU_header['TELESCOP'] == 'TAROT CALERN':
                                    telescope="TCA"
                                elif self.HDU_header['TELESCOP'] == 'TAROT CHILI':
                                    telescope="TCH"
                                elif self.HDU_header['TELESCOP'] == 'TAROT REUNION':
                                    telescope="TRE"
                                else:
                                    print("Error: export png to database\n")
                                    sys.exit()
                            try:
                                Insert_Image(db=telescope, ID=k['NUMBER'].item(), FoV=int(i),
                                         RA=k['X_WORLD'].item(), DE=k['Y_WORLD'].item(),
                                         fits=name_subfit_web, dt='gwevent', png=name_PNG_web)
                            except:
                                pass
                            
                            #write png to list in /tmp/
                            with open(path_PNG_list, 'a') as indexsite:
                                indexsite.write(name_PNG_web)
                                indexsite.write("\n")
                            """Export PNG for website without reference star and correct ra and dec """
                            # Input in X, Y pixel and polot directly but use projection wcs
                            if not os.path.exists(save_PNG_pathfile):
#                                Plot_Image_Candi(file_new, k["X_IMAGE"], k["Y_IMAGE"], save_PNG_pathfile)
                                Plot_Image_Candi_Solar_Obj(file_new, k["X_WORLD"].item(), k["Y_WORLD"].item(), save_PNG_pathfile)
                                
                                """Exprot subfits in case DSS field"""
                                dss_ima = LOAD_DSS_IMAGE()
                                dss_ima.Read_Log()
                                if dss_ima.link_ok:
                                    dss_ima.Load_Fits_DSS(save_subfit, float(k["X_WORLD"]), float(k["Y_WORLD"]))
                                else:
                                    fit= SUB_FITS(file_new)
                                    fit.Sub_Fits_2(k['X_WORLD'].item(), k['Y_WORLD'].item(), str(k['NUMBER'].item()), save_subfit)
                            else:
                                print("PNG exist: %s\n" %save_PNG_pathfile)
                        else:
                            print("No NEW to export PNG\n")
                            pass
                else:
                    pass

if __name__ == '__main__':
# Run test FoV
    home = os.path.expanduser('~')
    folder = os.path.join(home, "Documents/Data/GrandMa_alert")
    run = TAROT_FoV(folder)
    run.FoV_Table()
#    print("Max FoV: %s" %run.max_FoV)
#    print("List FoV\n", run.list_FoV)
    run.Predict_FoV()
#    print(run.table)
    run.ML_Candi()
#    print(run.table)
    run.Export_FoV_Candi()
    run.Export_PNG_Candi()
    print(run.table)