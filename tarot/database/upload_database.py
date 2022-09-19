#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 10:57:17 2020

@author: kanthanakorn
"""

import os
import shutil
import subprocess

def Copy_Folder(src, dest):
    # Copy folder
    runcode = ("cp -r %s %s" %(src, dest))
    run = subprocess.Popen(runcode, shell=True)
    out, err = run.communicate()
    print("error: %s" %err)

def Copy_File(src, dest):
    runcode = ("cp %s %s" %(src, dest))
    run = subprocess.Popen(runcode, shell=True)
    out, err = run.communicate()
    print("error: %s" %err)

def Move_File(src, dest):
    runcode = ("mv %s %s" %(src, dest))
    run = subprocess.Popen(runcode, shell=True)
    out, err = run.communicate()
    print("error: %s" %err)
    
def Move_Image_Folder(tele):
    #Image path for tca, tch, tch
    image_tca = '/opt/lampp/htdocs/tarot-tca/images'
    image_tch = '/opt/lampp/htdocs/tarot-tch/images'
    image_tre = '/opt/lampp/htdocs/tarot-tre/images'
    link_to_copy = '/tmp/keep_PNG/folder_PNG.txt'
    
    if tele == "TCA":
        if os.path.exists(image_tca):
            shutil.rmtree(image_tca)
            
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tca)
        else:
            pass
            
        with open(link_to_copy, 'r') as rfolder:
            path_data = rfolder.readlines()
        for i in path_data:
            Copy_Folder(i.rstrip(), image_tca)
            
    elif tele == "TCH":
        if os.path.exists(image_tch):
            shutil.rmtree(image_tch)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tch)
        else:
            pass
            
        with open(link_to_copy, 'r') as rfolder:
            path_data = rfolder.readlines()
        for i in path_data:
            Copy_Folder(i.rstrip(), image_tch)
            
    elif tele == "TRE":
        if os.path.exists(image_tre):
            shutil.rmtree(image_tre)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tre)
        else:
            pass
        
        with open(link_to_copy, 'r') as rfolder:
            path_data = rfolder.readlines()
        for i in path_data:
            Copy_Folder(i.rstrip(), image_tre)
    else:
        print("Error copy image to website TAROT\n")
        
def Drop_FOV_Folder(telescope):
    #Image path for tca, tch, tch
    image_tca = '/opt/lampp/htdocs/tarot-tca/images'
    image_tch = '/opt/lampp/htdocs/tarot-tch/images'
    image_tre = '/opt/lampp/htdocs/tarot-tre/images'
    
    if telescope == "TCA":
        if os.path.exists(image_tca):
            shutil.rmtree(image_tca)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tca)
        else:
            pass
    elif telescope == "TCH":
        if os.path.exists(image_tch):
            shutil.rmtree(image_tch)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tch)
        else:
            pass
    elif telescope == "TRE":
        if os.path.exists(image_tre):
            shutil.rmtree(image_tre)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tre)
        else:
            pass

def Drop_Galaxy_Folder(telescope):
    #Image path for tca, tch, tch
    image_tca = '/opt/lampp/htdocs/tarot-tca/galaxy'
    image_tch = '/opt/lampp/htdocs/tarot-tch/galaxy'
    image_tre = '/opt/lampp/htdocs/tarot-tre/galaxy'
    
    if telescope == "TCA":
        if os.path.exists(image_tca):
            shutil.rmtree(image_tca)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tca)
        else:
            pass
    elif telescope == "TCH":
        if os.path.exists(image_tch):
            shutil.rmtree(image_tch)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tch)
        else:
            pass
    elif telescope == "TRE":
        if os.path.exists(image_tre):
            shutil.rmtree(image_tre)
            #build again, other it will copy folder 1 to images
            os.makedirs(image_tre)
        else:
            pass

    
def Copy_File_To_Database(file_to_copy, fov_folder, telescope):
    if not isinstance(fov_folder, str):
        fov_folder = str(fov_folder)
    else:
        pass
   
    #Image path for tca, tch, tch
    image_tca = '/opt/lampp/htdocs/tarot-tca/images'
    image_tch = '/opt/lampp/htdocs/tarot-tch/images'
    image_tre = '/opt/lampp/htdocs/tarot-tre/images'
    
    if os.path.exists(file_to_copy):    
        if telescope == "TCA":
            folder = os.path.join(image_tca, fov_folder)
            if os.path.exists(folder):
                Copy_File(file_to_copy, folder)
            else:
                os.makedirs(folder)
                Copy_File(file_to_copy, folder)
                
        elif telescope == "TCH":
            folder = os.path.join(image_tch, fov_folder)
            if os.path.exists(folder):
                Copy_File(file_to_copy, folder)
            else:
                os.makedirs(folder)
                Copy_File(file_to_copy, folder)
                
        elif telescope == "TRE":
            folder = os.path.join(image_tre, fov_folder)
            if os.path.exists(folder):
                Copy_File(file_to_copy, folder)
            else:
                os.makedirs(folder)
                Copy_File(file_to_copy, folder)
        else:
            print("Error copy image to website TAROT\n")
    else:
        print("Not exists: %s" %file_to_copy)
