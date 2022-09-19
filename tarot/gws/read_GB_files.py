#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 11:33:30 2019

@author: noysena

Read information from xml to provide Burst time, delay time, Astrophysical
"""

import os
import json
import glob
#import copy
import shutil
import numpy as np
from datetime import datetime
import xml.etree.ElementTree as ET
from astropy.utils.data import download_file
from ligo.skymap.tool.ligo_skymap_plot import main as main1
from ligo.skymap.tool.ligo_area import main as main2
from ligo.skymap.tool.ligo_skymap_plot_volume import main as main3
from ligo.skymap.tool.ligo_distance import main as main4

#link = 'S190930s-1-Preliminary.xml'
#filename = "bayestar.fits.gz"
#read xml file to get event and astrophysics
def Read_XML(link):
    """ Computer burst time and preliminary recode time"""
    tree = ET.parse(link)
    root = tree.getroot()
    
    #Creat file (to see delay time)    
    for elem in root.findall("./Who"):
    #    print(i.tag, i.text)
         date = elem.find('Date').text
#         print("Creat Data : %s" %date)
         
         try:
             D1 = datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
         except(ValueError):
             D1 = datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f")
         except(AttributeError):
             dt1, _, us= D1.partition(".")
             D1 = datetime.strptime(dt1, "%Y-%m-%dT%H:%M:%S")
        
    
    #print(D1)
    
    #Burst Event    
    for elem in root.iter('ISOTime'):
        detected = elem.text
        dt0, _, us= detected.partition(".")
        event_date_time = datetime.strptime(dt0, "%Y-%m-%dT%H:%M:%S")
#        print("Burst event: %s" %detected)
        try:
            D0 = datetime.strptime(detected, "%Y-%m-%dT%H:%M:%S.%f")
        except(ValueError):
            dt0, _, us= detected.partition(".")
            D0 = datetime.strptime(dt0, "%Y-%m-%dT%H:%M:%S")
            
    #print(D0)         
    
    #Delay time
    delay = np.round(((D1 - D0).seconds /60.), 2)
#    print("Delay (min): %0.2f" %delay)
    
    #GW classification take from "p_astro.json link" which is more update
    info=[]
    for elem in root.iter('Param'):
        name = elem.get('name')
        value= elem.get('value')
    #    print(name, value)
        if name in ['GraceID', 'BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']:
            try:
                value = round(float(value)*100)
            except(ValueError):
                pass
            info.append(value)
#            print(name, value)
    return event_date_time, delay, info
        
# =============================================================================
# Skymap to get sq.deg for 50 and 90 contour levels
# =============================================================================
#package ligo.skymap file for data
#/home/noysena/.local/lib/python3.6/site-packages/ligo/skymap/tool/
#ligo_skymap_plot.py line : 155 - 156
#ligo_skymap_plot_volume.py line : 248+249

def Read_Fits_Gz(filename):
    
    name = os.path.basename(filename)
    
    source = name.replace('.gz', '.pdf')
    destination = os.path.dirname(filename)

    PLOT = False
    if PLOT:
        #only sq.deg with skymaps plot
        savefile1 = ('-o %s' %source)
        main1([filename, savefile1, '--annotate', '--contour', '50', '90'])
    
        #there is space in font of file
        fix_source = ' ' + source
        shutil.copy2(fix_source, destination)
        
        #distance with volume plot
        savefile2 = ('-o %s' %('volume_'+source))
        main3([filename, savefile2, '--annotate'])
        shutil.copy2(fix_source, destination)
        
        prob, dist = [[[50, 90], [0, 0]], [0, 0]]
        
    else:
        #only sq.deg without skymaps plot
        savefile1 = ('-o %s' %source)
        prob = main2([filename, savefile1, '--annotate', '--contour', '50', '90'])
#        print(prob[0], prob[1])
        
        #distance without volume plot
        savefile2 = ('-o %s' %('volume_'+source))
        dist = main4([filename, savefile2, '--annotate'])
#        print(dist[0], dist[1])
        
    return prob, dist

# =============================================================================
# Astrophysics
# =============================================================================
def Read_Astro(filename):
    with open(filename) as rj:
        rjson = json.load(rj)
    BSN = round(float(rjson["BNS"])*100)
    NSBH = round(float(rjson["NSBH"])*100)
    BBH = round(float(rjson["BBH"])*100)
    MassGap = round(float(rjson["MassGap"])*100)
    Terrestrial = round(float(rjson["Terrestrial"])*100)
    return BSN, NSBH, BBH, MassGap, Terrestrial
    
# =============================================================================
# search specific file
# =============================================================================
def Filter_Text(rjson, text):
    ########################################################        
    k=[]
    for i, j in enumerate(rjson):
        if text in j:
            k.append(j)
        else:
            pass
    if k != []:
        file = rjson[k[0]]
    else:
        file = None
    return file

def Last_File(rjson, text):
    ########################################################        
    k=[]
    for i, j in enumerate(rjson):
        if text in j:
            k.append(j)
        else:
            pass
    if k != []:
        file = rjson[k[-1]]
    else:
        file = None
#    print(rjson[k[0]])
    return file
 
def Check_File(file):
    try:
        is_file = glob.glob(file)[0]
    except(IndexError):
        is_file = "NA"
        print("No: %s" %file)
    return is_file
# =============================================================================
# Download file    
# =============================================================================
def Preliminary_Load(superevent, rjson):
    #xml 1-preliminary, the first list in json
    try:
        xml_name = superevent + "-1-Preliminary.xml"
        xml_link = rjson[xml_name]
        if os.path.exists(xml_name):
            print("File exists: %s" %xml_name)
        else:
            file_xml = download_file(xml_link, cache=True)
            shutil.copy2(file_xml, xml_name)
    except(KeyError):
        print("No Preliminary.xml")

def Initial_Load(superevent, rjson):
    xml_link = Filter_Text(rjson, "Initial.xml")
    if xml_link != None:
        xml_name = os.path.basename(xml_link)
        if os.path.exists(xml_name):
            print("File exists: %s" %xml_name)
        else:
            file_xml = download_file(xml_link, cache=True)
            shutil.copy2(file_xml, xml_name)
    else:
        print("No initial.xml")

def Update_Load(superevent, rjson):
    xml_link = Filter_Text(rjson, "Update.xml")
    if xml_link != None:
        xml_name = os.path.basename(xml_link)
        if os.path.exists(xml_name):
            print("File exists: %s" %xml_name)
        else:
            file_xml = download_file(xml_link, cache=True)
            shutil.copy2(file_xml, xml_name)
    else:
        print("No update.xml")

def Json_Load(superevent, rjson):
    try:
        json_link = rjson["p_astro.json"]
        file_name_json = superevent + "-p_astro.json"
        if os.path.exists(file_name_json):
            print("File exists: %s" %file_name_json)
        else:
            file_json = download_file(json_link, cache=True)
            shutil.copy2(file_json, file_name_json)
    except(KeyError):
        print("No -p_astro.json")
        
def Bayerstar_Load(superevent, rjson):
    try:
        fits_bayestar_link = rjson["bayestar.fits.gz"]
        file_name_fits = superevent + "-bayestar.fits.gz"
        if os.path.exists(file_name_fits):
            print("File exists: %s" %file_name_fits)
        else:
            file_fits_bayestar_load = download_file(fits_bayestar_link, cache=True)
            shutil.copy2(file_fits_bayestar_load, file_name_fits)
    except(KeyError):
        print("No bayerstar")

def Bayerstar_Update_Load(superevent, rjson):
    bay_update_link = Last_File(rjson, "bayestar.fits.gz")
    if bay_update_link != None:
        upif_load_name = os.path.basename(bay_update_link)
        fits_name = upif_load_name[:-2] #upif_load_name.replace("fits.gz,1", "fits.gz")
    
        upfi_name = superevent + "-update-" + fits_name
        if os.path.exists(upfi_name):
            print("File exists: %s" %upfi_name)
        else:
            file_upfi = download_file(bay_update_link, cache=True)
            shutil.copy2(file_upfi, upfi_name)
    else:
        print("No update Bayerstar")

def Lal_Load(superevent, rjson):
    lal_update_link = Last_File(rjson, "LALInference.fits.gz")
    if lal_update_link != None:
        upif_load_name = os.path.basename(lal_update_link)
        fits_name = upif_load_name[:-2]
        #the last update will be apply for lal, not the same as bayser
        upfi_name = superevent + "-" + fits_name
        if os.path.exists(upfi_name):
            print("File exists: %s" %upfi_name)
        else:
            file_fits_lal_load = download_file(lal_update_link, cache=True)
            shutil.copy2(file_fits_lal_load, upfi_name)
    else:
        print("No lal")
            
def Cwb_Load(superevent, rjson):
    cwb_link = Last_File(rjson, "cWB.fits.gz")
    if cwb_link != None:
        file_fits_cWB = superevent + "-cWB.fits.gz"
        if os.path.exists(file_fits_cWB):
            print("File exists: %s" %file_fits_cWB)
        else:
            file_fits_cWb_load = download_file(cwb_link, cache=True)
            shutil.copy2(file_fits_cWb_load, file_fits_cWB)
    else:
        print("No cWB") 
             
        
def Download_File(superevent):
    superevent_json = superevent + ".json"
    link = "https://gracedb.ligo.org/api/superevents/" + superevent + "/files/?format=json"
    print(link)
    if os.path.exists(superevent_json):
        print("File Superevent json exists: %s" %superevent_json)
        filejson = superevent_json
    else:
        filejson = download_file(link)
        shutil.copy2(filejson, superevent_json)
    
    with open(filejson, 'r') as rj:
        rjson = json.load(rj)

    #load preliminary
    Preliminary_Load(superevent, rjson)
    #load initial
    Initial_Load(superevent, rjson)
    
    #load update
    Update_Load(superevent, rjson)
    
    #load json
    Json_Load(superevent, rjson)
    
    #load bayerstar
    Bayerstar_Load(superevent, rjson)
    
    #load bayerstar update
    Bayerstar_Update_Load(superevent, rjson)
    
    #load lal
    Lal_Load(superevent, rjson)
    
    #load cWB
    Cwb_Load(superevent, rjson)
    
# =============================================================================
# MAIN {id:[pre, ini, update, json, fits, fits_update]}
# =============================================================================
if __name__ == "__main__":
    dpath = '/home/noysena/Documents/Code/Code_O3/All_GWs/'
    with open('dict-gw.txt', 'r') as rj:
        superevent = json.load(rj)
    
#    o3 = open("O3_info.txt", 'w')
    
    #generate Event and burst time from file preliminary and initail
    if True:    
        for i, j in superevent.items():    
            url = os.path.join(dpath, superevent[i]['preliminary'])
            try:
                burst, delay, csf = Read_XML(url)
            except(FileNotFoundError, IsADirectoryError):
                url = os.path.join(dpath, superevent[i]['initial'])
                burst, delay, csf = Read_XML(url)
            print("Event %s %s" %(i, burst))               
#        
    
#    for i in superevent:
#        Download_File(i)
#        
#        pre = i + "*Preliminary.xml"
#        try:
#            pre_file = glob.glob(pre)[0]
#        except(IndexError):
#            print("No: %s" %pre)
#            pre_file = "NA"
#            pass
#            #No preliminary exist, so use -1-initial.xml instead
#            print("No Preliminary, use Initial instead")
#            pre_ini = i + "*Initial.xml"
#            pre_file = glob.glob(pre_ini)[0]
#        
#        
#        ini = i + "*Initial.xml"
#        ini_file = Check_File(ini)
#            
#        upd = i + "*Update.xml"
#        upd_file = Check_File(upd)
#        
#        jso = i + "*p_astro.json"
#        jso_file = Check_File(jso)
#        
#        bay = i + "*bayestar.fits.gz"
#        bay_file = Check_File(bay)
#        
#        lal = i + "*LALInference.fits.gz"
#        lal_file = Check_File(lal)
#
#        cwb = i + "*cWB.fits.gz"
#        cwb_file = Check_File(cwb)
#        
#        upfi= i + "*update-bayestar.fits.gz"
#        upfi_file = Check_File(upfi)
#        
#        if pre_file != "NA":
#            pre_detect, pre_delay, pre_info = Read_XML(pre_file)
#        else:
#            #Initial file is alwals produced
#            print("No Prelimitary, use initial instead")
#            pre_detect, pre_delay, pre_info = Read_XML(ini_file)
#
#        if upd_file != "NA":
#            upd_detect, upd_delay, upd_info = Read_XML(upd_file)
#        else:
#            print("No Update, Use initial")
#            upd_detect, upd_delay, upd_info = Read_XML(ini_file)
#        
#        if jso_file != "NA":
#            BSN, NSBH, BBH, MassGap, Terrestrial = Read_Astro(jso_file)
#        else:
#            BSN, NSBH, BBH, MassGap, Terrestrial = [0, 0, 0, 0, 0]
#            print("No astro-json %s: " %jso)
#        
#        if bay_file != "NA":
#            prob_bay, dist_bay = Read_Fits_Gz(bay_file)
#        else:
#            prob_bay, dist_bay = [[[50, 90], [0, 0]], [0, 0]]
#            print("No Bayerstar %s: " %bay_file)
#        
#        #pudate fits
#        if lal_file != "NA":
#            print("Update maps with LAL")
#            prob, dist = Read_Fits_Gz(lal_file)
##        elif cwb_file != "NA":
##            print("Update maps with cWB")
##            prob, dist = Read_Fits_Gz(cwb_file)
#        elif upfi_file != "NA":
#            print("Update maps with Bayerstar")
#            prob, dist = Read_Fits_Gz(upfi_file)
#        else:
#            print("No LAL & Update-Baystar")
#            prob, dist = [[[50, 90], [0, 0]], [0, 0]]
#   
#        print(pre_info[0], pre_detect, pre_delay,
#            prob_bay[1][0], prob_bay[1][1], dist_bay[0], dist_bay[1], upd_delay,
#            prob[1][0], prob[1][1], dist[0], dist[1],
#            BSN, NSBH, BBH, MassGap, Terrestrial)
#        
#        o3_data = ("%s & %s & %2.2f & %d & %d & $%d \pm %d$ & %d & %d & %d & $%d \pm %d$ & %d & %d & %d & %d & %d \\\\ \n" %(pre_info[0], pre_detect, pre_delay,
#                   prob_bay[1][0], prob_bay[1][1], dist_bay[0], dist_bay[1], upd_delay,
#                   prob[1][0], prob[1][1], dist[0], dist[1],
#                   BSN, NSBH, BBH, MassGap, Terrestrial))
#        o3.write(o3_data)
#        
#        print("+"*80)
#    o3.close()
