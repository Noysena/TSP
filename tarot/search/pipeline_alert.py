#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:33:28 2019

@author: kanthanakorn
"""

# =============================================================================
# SMS pipeline report
# =============================================================================
import requests
import os

"""Need file in /tmp/status.txt and write status down by Pipeline;
SMS_alert("/tmp/status.txt", "Hellow%20AstroNu%20TSP%20Pipeline")"""

def SMS_alert(log, keywords):
    if os.path.exists(log):
        with open(log, 'r') as status:
            read_status = status.read()
            if read_status.rstrip() == 'GWalert':
                user = "27094265"
                pwd = "jST8oTFvsJXVHs"
                text = ("https://smsapi.free-mobile.fr/sendmsg?user=%s&pass=%s&msg=%s" %(user, pwd, keywords))
                requests.post(text, data = {'key':'value'})
            else:
                pass
    else:
        print("No log: %s" %log)
        
def SMS_direct_alert(keywords):
    user = "27094265"
    pwd = "jST8oTFvsJXVHs"
    text = ("https://smsapi.free-mobile.fr/sendmsg?user=%s&pass=%s&msg=%s" %(user, pwd, keywords))
    requests.post(text, data = {'key':'value'})
