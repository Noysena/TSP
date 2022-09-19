#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 13:06:13 2019

@author: kanthanakorn

This funtion is to read link GrandMa alert export by CADOR.
"""
import os
import sys
import copy
import shutil
import requests
from bs4 import BeautifulSoup
from astropy.table import Table
from urllib.request import urlopen
from tarot.utility.test_run import Test_Run
from astropy.utils.data import download_file

class TAROT_Retrieve:
    def MAIN():
        return("Searching TAROT log for GW alert")
        
    def __init__(self):
        self.link_url = ''
        self.read_log = ''
        self.alert_path = "/mnt/tarot/onsite_storage/al_onsite"
        self.alert_path_TCA = "/mnt/tarot/onsite_storage/al_onsite/TCA"
        self.alert_path_TCH = "/mnt/tarot/onsite_storage/al_onsite/TCH"
        self.alert_path_TRE = "/mnt/tarot/onsite_storage/al_onsite/TRE"
        self.main_link_tca = "http://tca4.tarotnet.org/ros/alert_grandma/"
        self.main_link_tch = "http://tch4.tarotnet.org/ros/alert_grandma/"
        self.main_link_tre = "http://lesmakes.dlinkddns.com:8084/ros/alert_grandma/"
        self.fits_log = "/tmp/grandma_alert_files.log"
        self.GrandMaAlert = False
        self.test = False
    
    def Read_Event_Page(self, link):
        self.gws = []
        open_link = urlopen(link)
        soup = BeautifulSoup(open_link, 'html.parser')
        for i in soup.find_all('a'):
            event_link = i.get('href')
            if event_link[0] == 'S':
                self.gws.append(event_link)
            else:
                pass
        
                
    def Check_Log_link(self, tele, event_html):
        self.test = Test_Run()
        print("Test status %s\nRun on machine: %s" %(self.test, os.getcwd()))
        
        if self.test:
            self.alert_path = "/tmp"
        else:
            pass
        
        # event_html = str(gwevent) + '.html'
        
        if tele == "TCA":
            self.link_url = self.main_link_tca + event_html
        elif tele == "TCH":
            self.link_url = self.main_link_tch + event_html
        elif tele == "TRE":
            self.link_url = self.main_link_tre + event_html
        else:
            self.link_url = False
            
        # return self.link_url
        
    def Read_Log(self, tele):
        """ Read log on observing night"""
        try:
            Check_link = requests.get(self.link_url, timeout=60)
        except:
            Check_link = False
            print("%s --> Error to get link " %self.link_url)
            pass

        if Check_link != False:
            Log_today = urlopen(self.link_url)
            print ("\nLink: %s --> %s\n" %(self.link_url, Check_link.reason))
            self.read_log = copy.copy(Log_today.read())
            
            # Read href link and download
            self.soup = BeautifulSoup(self.read_log, 'html.parser')
            
            try:
                table = self.soup.body.pre
            except:
                table = self.soup.find('pre')
                
            table_text = copy.copy(table)
            
            if tele == "TCA":
                list_fits = "/tmp/TCA_FITS_files.txt"
            elif tele == "TCH":
                list_fits = "/tmp/TCH_FITS_files.txt"
            elif tele == "TRE":
                list_fits = "/tmp/TRE_FITS_files.txt"
            else:
                list_fits = "/tmp/Tarot_FITS_files.txt"
                
            write_list_fits = open(list_fits, 'w')
            
            for i in table_text:
                try:
                    write_list_fits.write(i.contents)
                except(AttributeError):
                    write_list_fits.write(i)
                except(TypeError):
                    write_list_fits.write(i.contents[0])
            
            write_list_fits.close()
            
            self.table = Table.read(list_fits, format='ascii')
            self.table.write("/tmp/Table_Lits_GW_FITS.dat", format='csv',
                             overwrite=True)
            
        # return self.table
    
    def Search_for_Fits(self):
        self.keep = []
        for link in self.soup.find_all('a'):
            self.keep.append(link.get('href'))
            print(link.get('href'))
            
    def Write_File_Link_to_Local(self, tele, event, file_link):
        
        if tele == "TCA":
            save_img = os.path.join(self.alert_path_TCA, event, os.path.basename(file_link))
        elif tele == "TCH":
            save_img = os.path.join(self.alert_path_TCH, event, os.path.basename(file_link))
        elif tele == "TRE":
            save_img = os.path.join(self.alert_path_TRE, event, os.path.basename(file_link))
        else:
            save_img = os.path.join(self.alert_path, event, os.path.basename(file_link))
            
        if os.path.exists(os.path.dirname(save_img)):
            print("Path exist: %s" %os.path.dirname(save_img))
        else:
            try:
                os.makedirs(os.path.dirname(save_img))
            except(FileExistsError):
                print("Path exist: %s" %os.path.dirname(save_img))
        
        if os.path.exists(save_img):
            print("File exist... done: %s" %os.path.basename(save_img))
        else:
            loadfile = download_file(file_link)
            shutil.copy2(loadfile, save_img)
            print(save_img)
            
    def Load_Fits(self, tele, gwevent):
        
        """ Load FITS to folder """
        for i, j in enumerate(self.soup.find_all('a')):
            if self.table[i]['WCSpb'] == 1:
                pass
            elif self.table[i]['sname'][0:9] != gwevent:
                pass
            else:
#                print("Load %s" %self.table[i])
                try:
                    TAROT_Retrieve.Write_File_Link_to_Local(self,
                                                            tele,
                                                            self.table[i]['sname'][0:9],
                                                            j.get('href'))
                except:
                    e = sys.exc_info()
                    print(e)


if __name__ == '__main__':
    if False:
        TCA = TAROT_Retrieve()
        TCA.Read_Event_Page(TCA.main_link_tca)
    
        for link_gw in TCA.gws:
            if len(link_gw) == 13:
                gw = link_gw.replace('.html', '_')
            else:
                gw = link_gw.replace('.html', '')
            
        TCA.Check_Log_link('TCA', link_gw)
        TCA.Read_Log('TCA')
        TCA.table

    if True:
        TCH = TAROT_Retrieve()
        TCH.Read_Event_Page(TCH.main_link_tch)

        for link_gw in TCH.gws:
            if len(link_gw) == 13:
                gw = link_gw.replace('.html', '_')
            else:
                gw = link_gw.replace('.html', '')

            TCH.Check_Log_link('TCH', link_gw)
            TCH.Read_Log('TCH')
            TCH.table
