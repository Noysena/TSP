#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# SMS pipeline gw alert
# =============================================================================
import os
import sys
import requests
#from datetime import datetime

# =============================================================================
# Free SMS alert
# =============================================================================
def Go_Free(keywords):
    user = "27094265"
    pwd = "jST8oTFvsJXVHs"                                                                                                                                  
    text = ("https://smsapi.free-mobile.fr/sendmsg?user=%s&pass=%s&msg=%s" %(user, pwd, keywords))                                                                               
    requests.post(text, data = {'key':'value'})

def SMS_direct_alert(keywords):
    if os.path.exists("sent.log"):
        with open("sent.log" , 'r') as slog:
            rsent = slog.read()
        if keywords in rsent:
            pass
        else:
            Go_Free(keywords)
            with open("sent.log", 'a') as slog:
                slog.write("%s\n" %keywords)
    else:
        with open("sent.log" , 'w') as slog:                                                                                                                                             
            slog.write("%s\n" %keywords)
        Go_Free(keywords)


# =============================================================================
# Read from graceGB
# =============================================================================
def Read_Events(link_url):
    """ Read log on observing night"""
    try:
        get_link = requests.get(link_url, timeout=3)
    except:
        get_link = None

    if get_link != None and get_link.ok:
        retrieve_json = get_link.json()
        #retrieve_json['numRows'] #Will tell that there is update event
    else:
        print("%s --> Error to convert link to JSON: " %link_url)
        sys.exit()

    superevent_all = retrieve_json['superevents']
    superevent_links = retrieve_json['links']
    
    try:
        last_event = superevent_all[0]['superevent_id']
        last_event_time = superevent_all[0]['created']
    except:
        last_event = ''

    return superevent_all, superevent_links, last_event, last_event_time

def gwlog(gw):
    recod_gw_event = "/tmp/gwevents.txt"
    if os.path.exists(recod_gw_event):
        with open(recod_gw_event, 'r') as gwtxt:
            readgw = gwtxt.read()
        if gw not in readgw:
            SMS_direct_alert("%s" %gw)
        else:
            pass

        if gw not in readgw:
            with open(recod_gw_event, 'a') as gwtxt:
                gwtxt.write("%s\n" %gw)
        else:
            pass
    else:
        with open(recod_gw_event, 'w') as gwtxt:
            gwtxt.write("%s\n" %gw)

        SMS_direct_alert("%s" %gw)

def True_Event(all_events):
    #GraceGB puts Superevent to [0] but we think that will keep read all S in a day
    for event in all_events:
        if event['superevent_id'][0:3] == 'S19':
            gw_event = event['superevent_id']
            gw_event_time = event['created']
            gwlog(gw_event)
            return gw_event, gw_event_time
        else:
            gw_event = None
            gw_event_time = None
    return gw_event, gw_event_time

def Test_Event(all_events, ms):
    for event in all_events:
        if event['superevent_id'] == ms:
            gw_event = event['superevent_id']
            gw_event_time = event['created']
            gwlog(gw_event)
            return gw_event, gw_event_time
        else:
            gw_event = None
            gw_event_time = None
    return gw_event, gw_event_time
        

if __name__ == '__main__':
#Link gracegb
    link = "https://gracedb.ligo.org/apiweb/superevents/?format=json"
#test = SMS_direct_alert("S190722a")

#Set iso and test event?
    iso_time = True
    test_event = False
    
    superevent_all, superevent_links, last_event, last_event_time = Read_Events(link)

    if test_event:
        ms = 'MS191226r'
        gw_event, gw_event_time = Test_Event(superevent_all, ms)
    else:
        gw_event, gw_event_time = True_Event(superevent_all)
