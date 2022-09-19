#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 10:06:38 2019

@author: noysena
"""

import os
import json
import shutil
import numpy as np
from datetime import datetime
import xml.etree.ElementTree as ET
from astropy.utils.data import download_file
from ligo.skymap.tool.ligo_skymap_plot import main as main1
from ligo.skymap.tool.ligo_area import main as main2
from ligo.skymap.tool.ligo_skymap_plot_volume import main as main3
from ligo.skymap.tool.ligo_distance import main as main4

# =============================================================================
# Function to check and download file from superevent
# =============================================================================
#When we have update file, this will help to download the last update
def Filter_Text(rjson, text):
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
        
    return xml_name

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
    return xml_name

#Download file        
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
    preli_name = Preliminary_Load(superevent, rjson)
    #load initial
    initi_name = Initial_Load(superevent, rjson)
    
    return preli_name, initi_name

#read xml file to get event and astrophysics    
def Read_XML(link):
    """ Computer burst time and preliminary recode time"""
    tree = ET.parse(link)
    root = tree.getroot()
    
    #Creat file (to see delay time)    
    for elem in root.findall("./Who"):
         date = elem.find('Date').text
         
         try:
             D1 = datetime.strptime(date, "%Y-%m-%dT%H:%M:%S")
         except(ValueError):
             D1 = datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f")
         except(AttributeError):
             dt1, _, us= D1.partition(".")
             D1 = datetime.strptime(dt1, "%Y-%m-%dT%H:%M:%S")
    #Burst Event    
    for elem in root.iter('ISOTime'):
        #convert date time in other 'astropy.Time('2019-04-08T18:18:02')'
        #import parser & dateutil.parser.parse('2019-04-08T18:18:02')
        detected = elem.text
        dt0, _, us= detected.partition(".")
        event_time = datetime.strptime(dt0, "%Y-%m-%dT%H:%M:%S")
#        print("Burst Event: %s" %event_time.isoformat())
#        print(event_time.isoformat())
        try:
            D0 = datetime.strptime(detected, "%Y-%m-%dT%H:%M:%S.%f")
        except(ValueError):
            dt0, _, us= detected.partition(".")
            D0 = datetime.strptime(dt0, "%Y-%m-%dT%H:%M:%S")
            
    #Delay time in minute
    delay = np.round(((D1 - D0).seconds /60.), 2)
    
    #GW classification take from "p_astro.json link" which is more update
    classification=[]
    for elem in root.iter('Param'):
        name = elem.get('name')
        value= elem.get('value')
        if name in ['GraceID', 'BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']:
            try:
                value = round(float(value)*100)
            except(ValueError):
                pass
            classification.append(value)
    return event_time, delay, classification

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
        prob = main2([filename, '--annotate', '--contour', '50', '90'])
        
        #distance without volume plot
        dist = main4([filename, '--annotate'])
        
    return prob, dist

if __name__ == '__main__':
# =============================================================================
#     Load current event and upload to database by Insert_Gws
# =============================================================================
    if True:
        from tarot.database.TSP_SQL import Insert_Gws
        import sys
#        last_gw = 'S191213g'
        last_gw = sys.argv[1]
        preli_name, init_name = Download_File(last_gw)
        print(preli_name, init_name)
        
        if os.path.isfile(preli_name):
            event, delay, csf = Read_XML(preli_name)
        else:
            event, delay, csf = Read_XML(init_name)
            
        superevent = [csf[0], event.isoformat().replace("T", " "), str(csf[1:6]), 0, 0, 0]
        print(superevent)
        Insert_Gws(superevent)
        
    else:
        pass
    
# =============================================================================
# Get Burst time event from Gws.txt (all event)
# =============================================================================
    if False:
        from tarot.gws.gws_list import Gws_Event
        gwfolder, gwname = Gws_Event()
        gws = {}
        for gw in gwname:
            preli_name, init_name = Download_File(gw)
            print(preli_name, init_name)
            
            if os.path.isfile(preli_name):
                event, delay, csf = Read_XML(preli_name)
            else:
                event, delay, csf = Read_XML(init_name)
                
            gws[gw] = event.isoformat()
            
        print(gws)
    else:
        pass
    
# =============================================================================
# Get Brust time event form specific event
# =============================================================================
    if False:
        import sys
        gw = sys.argv[1]
        preli_name, init_name = Download_File(gw)
        print(preli_name, init_name)
        
        if os.path.isfile(preli_name):
            event, delay, csf = Read_XML(preli_name)
        else:
            event, delay, csf = Read_XML(init_name)
            
        print(gw, event, delay, csf)
    else:
        pass
# =============================================================================
#     Load dict form /tmp/gws_burst_time.txt and export for event.php
#    if file is not in /tmp/it will laod form tarot/gws/
# =============================================================================
    if False:
        #dict-gw.txt generate form theme dict_gw.txt, therefore it must be updated.
        with open('dict-gw.txt') as gwsj:
            superevent = json.load(gwsj)
        
        dpath = '/home/noysena/Documents/Code/Code_O3/All_GWs/'
        gws_dict = {}
        if True:
            for i, j in superevent.items():
                url = os.path.join(dpath, superevent[i]['preliminary'])
                try:
                    event, delay, csf = Read_XML(url)
                except(FileNotFoundError, IsADirectoryError):
#                    print("No file")
                    url = os.path.join(dpath, superevent[i]['initial'])
                    event, delay, csf = Read_XML(url)
                #buil dictionaries
                gws_dict[i] = {'burst' : event.isoformat()}
                
                #update dictinaries
                gws_dict[i].update({'delay' : delay})
                gws_dict[i].update({'classification' : csf})
        dumps_gws = json.dumps(gws_dict)
        
        with open('gws_burst_time.txt', 'w') as rgws:
            rgws.write(dumps_gws)
        print(gws_dict)
        
        for i, j in gws_dict.items():
            
            print("('%s', '%s', '%s', 0, 0, 0)," %(i, j['burst'].replace("T"," "), j['classification'][1:6]))
#            insert_gws_table = ("('%s', '%s', '%s', 0, 0, 0)," %(i, j['burst'].replace("T"," "), j['classification'][1:6]))
#            Insert_Gws(insert_gws_table)
    else:
        pass
    
    #export file
    if False:
        with open('gws_burst_time.txt') as bt:
            burst_time = json.load(bt)
        print(burst_time)
    else:
        pass
    
