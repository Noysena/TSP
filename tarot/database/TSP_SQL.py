#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:18:33 2019

@author: noysena

Write database from TSP
"""
import datetime
import subprocess
import mysql.connector
from mysql.connector import errorcode

# =============================================================================
# Connect to exist DB
# =============================================================================
def Call_Database(tb=''):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %tb)
        
        mycursor = mydb.cursor()
        
        if True:
            mycursor.execute("SHOW DATABASES")
            
            print("Databases")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
        
    mydb.close() 

# =============================================================================
# Call TABLE
# =============================================================================
def Call_Table(db='Pydatabase'):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
        if True:
            mycursor.execute("SHOW TABLES")
            
            print("Databases")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
        
    mydb.close()
    
# =============================================================================
# Describes table
# =============================================================================
def Describe_Table(db='Pydatabase', tb='GWs'):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
        if True:
            mycursor.execute("DESCRIBE %s" %tb)
            
            print("Databases")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
        
    mydb.close()
# =============================================================================
# Build database
# =============================================================================
def Create_Database(nameit):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!"
        )
        
        mycursor = mydb.cursor()
        
        if False:
            mycursor.execute("SHOW DATABASES")
            
            print("Databases")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    mycursor.execute("CREATE DATABASE %s" %nameit)
    mydb.commit()
    mydb.close()
 
    ## =============================================================================
## Create table superevent --> index.php
## =============================================================================
def Create_Table(tb, db = "Pydatabase"):
    """
    GW events were created and if it is exists, it will pass.
    """
    print("Create TABLE: %s" %tb)
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    if False:
        mycursor.execute("SHOW TABLES")
        print("TABLES")
    
        tables = []
        for i, j in enumerate(mycursor):
            print(i, j)
            tables.append(j[0])
    
    columns = "ID int NOT NULL UNIQUE, FoV int not null, png VARCHAR(255)DEFAULT NULL,\
    gif VARCHAR(255) DEFAULT NULL, lightcurve VARCHAR(255) DEFAULT NULL, \
    data VARCHAR(255) DEFAULT NULL, RA float(6,3) DEFAULT NULL, \
    DE float(6,3) DEFAULT NULL, fits VARCHAR(255) DEFAULT NULL"
    mycursor.execute("CREATE TABLE %s (%s)" %(tb, columns))
    mydb.commit()
    mydb.close()

# =============================================================================
# Create table for GW list --> event.php
# =============================================================================
def Create_Gws_Table(tb="GWs"):
    print("Create TABLE: %s" %tb)
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    if False:
        mycursor.execute("SHOW TABLES")
        print("TABLES")
    
        tables = []
        for i, j in enumerate(mycursor):
            print(i, j)
            tables.append(j[0])
        
    columns = "`ID` int not null AUTO_INCREMENT,`event` varchar(10) DEFAULT NULL,\
    `burst` datetime DEFAULT NULL,`classification` varchar(20) DEFAULT NULL, \
    `tca` int(9) DEFAULT NULL, `tch` int(9) DEFAULT NULL,`tre` int(9) DEFAULT NULL,\
    PRIMARY KEY (`ID`)"
#    mycursor.execute("CREATE TABLE %s (%s) ENGINE=InnoDB DEFAULT CHARSET=latin1" %(event, columns))
    try:
        mycursor.execute("CREATE TABLE %s (%s)" %(tb, columns))
        mydb.commit()
    except:
       print("Error or Table exists: %s" %tb)
    mydb.close()

# =============================================================================
# Create galaxy table    
# =============================================================================
def Create_Galaxy_Table(db="Pydatabase", tb="galaxy"):
    print("Create TABLE: %s on %s" %(tb, db))
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    if False:
        mycursor.execute("SHOW TABLES")
        print("TABLES")
    
        tables = []
        for i, j in enumerate(mycursor):
            print(i, j)
            tables.append(j[0])
        
    columns = "`ID` int not null AUTO_INCREMENT,`FITS` varchar(73) DEFAULT NULL,\
    `PNG` varchar(73) DEFAULT NULL,`RA` float(6,3) DEFAULT NULL, `DEC` float(6,3) \
    DEFAULT NULL, `GIF` varchar(73) DEFAULT NULL, PRIMARY KEY (`ID`)"
#    mycursor.execute("CREATE TABLE %s (%s) ENGINE=InnoDB DEFAULT CHARSET=latin1" %(event, columns))
    try:
        mycursor.execute("CREATE TABLE %s (%s)" %(tb, columns))
        mydb.commit()
    except:
       print("Error or Table exists: %s" %tb)
    mydb.close()
# =============================================================================
# Create Alert User talbe
# =============================================================================
def Create_Alert_User_Table(tb="AlertUser"):
    print("Create TABLE: %s" %tb)
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    if False:
        mycursor.execute("SHOW TABLES")
        print("TABLES")
    
        tables = []
        for i, j in enumerate(mycursor):
            print(i, j)
            tables.append(j[0])
        
    columns = "`Alert` varchar(255) DEFAULT NULL, `info` varchar(255) DEFAULT NULL"
#    mycursor.execute("CREATE TABLE %s (%s) ENGINE=InnoDB DEFAULT CHARSET=latin1" %(tb, columns))
#    mydb.commit()
    mycursor.execute("CREATE TABLE %s (%s)" %(tb, columns))
    mydb.commit()
    mydb.close()
    
# =============================================================================
# Drop database name
# =============================================================================
def Drop_Database(dbname):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!")
        
        mycursor = mydb.cursor()
        
        if False:
            mycursor.execute("SHOW DATABASES")
            
            print("Table")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    print("Drop Database %s" %dbname)
    mycursor.execute("DROP DATABASE %s" %dbname)
    mydb.commit()
    mydb.close()

# =============================================================================
# Drop tables
# =============================================================================
def Drop_Table(tb, db = "Pydatabase"):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
        if False:
            mycursor.execute("SHOW TABLES")
            print("Table")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    print("Drop table %s on Database %s" %(tb, db))
    mycursor.execute("DROP TABLE %s" %tb)
    mydb.commit()
    mydb.close()
    
# =============================================================================
# Change table name
# =============================================================================
def Alter_Table_Name(tb, new_tb, db = "Pydatabase"):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
        if False:
            mycursor.execute("SHOW TABLES")
            print("Table")
            for i, j in enumerate(mycursor):
                print("%s: %s" %(i,j))
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    print("Change table name %s to %s in Database %s" %(tb, new_tb, db))
    mycursor.execute("ALTER TABLE %s RENAME TO %s" %(tb, new_tb));
    mydb.commit()
    mydb.close() 
 
# =============================================================================
# Call GW evetn from function gw_event
# =============================================================================
#def Call_GW_Event():
#    if os.path.exists("/tmp/gwevents.txt"):
#        with open("/tmp/gwevents.txt", 'r') as rgw:
#            all_events = rgw.readlines()
#    else:
#        with open("gwevents.txt", 'r') as rgw:
#            all_events = rgw.readlines()
#    
#    event = str(all_events[-1])
#    return all_events, event

# =============================================================================
# kwargs function to assigne attribute, need it
# =============================================================================
def kwargs_fn(dic, dics):
    try:
        dic = dics[dic]
    except(KeyError):
        dic = None
    return dic

# =============================================================================
# Upload data to databease, when we have png, gif, lc, data.
# =============================================================================
def Insert_Candidate(db="TCA", **kwargs):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    ID = kwargs_fn('ID', kwargs)
    FoV = kwargs_fn('FoV', kwargs)
    png = kwargs_fn('png', kwargs)
    gif = kwargs_fn('gif', kwargs)
    lc = kwargs_fn('lc', kwargs)
    data = kwargs_fn('data', kwargs)
    RA = kwargs_fn('RA', kwargs)
    DE= kwargs_fn('DE', kwargs)
    fits = kwargs_fn('fits', kwargs)
    dt = kwargs_fn('dt', kwargs)
    
    for i, j in kwargs.items():
        print(i,j)
        
    sql = ("INSERT INTO `%s`(`ID`, `FoV`, `png`, `gif`, `lightcurve`, `data`, `RA`, `DE`, `fits`) VALUES ('%d', '%d', '%s', '%s', '%s', '%s', '%f', '%f', '%s')" %(dt, ID, FoV, png, gif, lc, data, RA, DE, fits))
    print(sql)
    mycursor.execute(sql)
    mydb.commit()
    mydb.close()

#Insert_Candidate(db="TCA", ID=737, FoV=1, png='images/png/1/IM_20190425_185709396_000001_86598208_1_737.PNG', gif='images/png/1/1_737.gif', lc='images/png/1/TLC_1_737.PNG', data='images/png/1/Photometry_1_737.dat', RA=244.453, DE=23.657, fits='images/png/1/DSS_1_737.fits', dt='gwevent') 
# =============================================================================
# Insert Candidate ID, FoV, png, RA, DEC, fits
# =============================================================================
def Insert_Image(db="TCA", **kwargs):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    ID = kwargs_fn('ID', kwargs)
    FoV = kwargs_fn('FoV', kwargs)
    png = kwargs_fn('png', kwargs)
    RA = kwargs_fn('RA', kwargs)
    DE= kwargs_fn('DE', kwargs)
    fits = kwargs_fn('fits', kwargs)
    dt = kwargs_fn('dt', kwargs)
    
    for i, j in kwargs.items():
        print(i,j)
        
    sql = ("INSERT INTO `%s`(`ID`, `FoV`, `png`, `RA`, `DE`, `fits`) VALUES ('%d', '%d', '%s', '%f', '%f', '%s')" %(dt, ID, FoV, png, RA, DE, fits))
    print(sql)
    mycursor.execute(sql)
    mydb.commit()
    mydb.close()
#Insert_Image(db="TCA", ID=237, FoV=2, RA=242.878, DE=23.597, fits='images/png/1/DSS_1_237.fits', dt='gwevent', png='images/png/1/IM_20190425_202556422_000001_86599608_1_237.PNG')
# =============================================================================
# Insert galaxy
# =============================================================================
def Insert_Galaxy_Image(db="TCA", dt="galaxy", **kwargs):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    fits = kwargs_fn('fits', kwargs)
    png = kwargs_fn('png', kwargs)
    RA = kwargs_fn('RA', kwargs)
    DE= kwargs_fn('DE', kwargs)
    gif = kwargs_fn('gif', kwargs)
    
    for i, j in kwargs.items():
        print(i,j)
        
    sql = ("INSERT INTO `galaxy`(`FITS`, `PNG`, `RA`, `DEC`, `GIF`) VALUES ('%s', '%s', '%f', '%f', '%s')" %(fits, png, RA, DE, gif))
    print(sql)
    mycursor.execute(sql)
    mydb.commit()
    mydb.close()
# =============================================================================
# update gif, lightcurve, data use this after Insert_Image(), refer ID and FoV
# =============================================================================
def Update_Photometry(db="TCA", **kwargs):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "%s" %db)
        
        mycursor = mydb.cursor()
        
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    
    ID = kwargs_fn('ID', kwargs)
    FoV = kwargs_fn('FoV', kwargs)
    gif = kwargs_fn('gif', kwargs)
    lc = kwargs_fn('lc', kwargs)
    data = kwargs_fn('data', kwargs)
#    RA = kwargs_fn('RA', kwargs)
#    DE= kwargs_fn('DE', kwargs)
    dt = kwargs_fn('dt', kwargs)
    
#    for i, j in kwargs.items():
#        print(i,j)
    
    sql = ("UPDATE `%s` SET `gif`='%s',`lightcurve`='%s',`data`='%s' WHERE ID='%d' and FoV='%d' " %(dt, gif, lc, data, ID, FoV))   
    
    #has RA and DEC but we need to round float
#    sql = ("UPDATE `%s` SET `gif`='%s',`lightcurve`='%s',`data`='%s' WHERE ID='%d' and FoV='%d' and RA='%f' and DE='%f'" %(dt, gif, lc, data, ID, FoV , RA, DE))
    
    print(sql)
    mycursor.execute(sql)
    mydb.commit()
    mydb.close()

#Update_Photometry(db="TCA", ID=237, FoV=2, dt='gwevent', gif='images/png/1/1_237.gif', lc='images/png/1/TLC_1_237.PNG', data='images/png/1/Photometry_1_237.dat')
# =============================================================================
# Inser Alert User
# =============================================================================
def Insert_Alert(alert_messege, superevent):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    sql = "INSERT INTO `AlertUser`(`Alert`, `Info`) VALUES ('%s','%s')" %(alert_messege, superevent)
    print(sql)
    mycursor.execute(sql)
    mydb.commit()
    
# =============================================================================
# update AlertUser        
# =============================================================================
def Update_Alert(alert_messege, superevent):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    sql = "UPDATE AlertUser SET Alert='%s', Info='%s' WHERE 1" %(alert_messege, superevent)
    print(sql)
    mycursor.execute(sql)
    mydb.commit()

# =============================================================================
# Insert table Gws --> event.php
# To upload all tiles use 'gws_info.py' and run, !change name of event!
    
# Or load  'gws_burst_time.txt' with json import.
#for i,j in gws_dict.items():
#    Insert_Gws([i, j['burst'].replace("T"," "), str(j['classification'][1:6]), 0, 0, 0])
# =============================================================================
def Insert_Gws(upl):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    sql = ("INSERT INTO `GWs`(`event`, `burst`, `classification`, `tca`, `tch`, `tre`) \
           VALUES ('%s', '%s', '%s', '%d', '%d', '%d')" %(upl[0], upl[1], upl[2], upl[3],
                                              upl[4], upl[5]))
    print(sql)
    mycursor.execute(sql)
    mydb.commit()

# =============================================================================
# # Update table Gws --> event.php
#    To update tiles, read list of dict from file TAROT_tiles.txt
#    which is generated by read_tiles.sh in 'TSP' folder
#    just copy text and take ',' in the last line out and past in terminal
#    Don't forget to get rid of "_"
#    copy list and pass then use loop as below
#    all_tiles = [{'S191205ah': [0, 12, 0, 'S191205ah']},{'S191204r_': [0, 206, 0, 'S191204r_']}]
#        for i in all_tiles:
#            for m, n in i.items():
#                print("Update tiles: %s" %m)
#                Update_Gws(n)    
# =============================================================================
def Update_Gws(tarot_tiles):
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    sql = ("UPDATE GWs SET tca=%d, tch=%d, tre=%d WHERE event='%s'" %(tarot_tiles[0],
                                                                    tarot_tiles[1],
                                                                    tarot_tiles[2],
                                                                    tarot_tiles[3]))
    print(sql)
    mycursor.execute(sql)
    mydb.commit() 
# =============================================================================
# Call event name and burst time for photometry
# =============================================================================
def Call_Event_Info():
    try:
        mydb = mysql.connector.connect(
        host = "localhost",
        user = "root",
        passwd = "deckanpolar2!",
        database = "Pydatabase")
        
        mycursor = mydb.cursor()
            
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
            
    sql = ("SELECT event, burst FROM `GWs`")
    print(sql)
    mycursor.execute(sql)
    myresult = mycursor.fetchall()
    
    gwlist = [x for x in myresult]
    gw_dict = {}
    for i, j in gwlist:
        gw_dict[i] = j.isoformat()
        
    mydb.commit()
    
    return gw_dict

# =============================================================================
# Backup datebase
# =============================================================================
def Backup_Database(db="Pydatabase"):
    main_databases = ["Pydatabase", "TCA", "TCH", "TRE"]
    save_date = datetime.date.today().isoformat()
    print("What do you want to back up?")
    Call_Database()
    
    if db == 'All':
        for d in main_databases:
            savesql = d + "_" + save_date + ".sql"
            dumpsql = ("mysqldump --databases %s > %s" %(d, savesql))
            run = subprocess.Popen(dumpsql, stdout=subprocess.PIPE, shell=True)
            out, err = run.communicate()
            print(savesql)
    else:
        savesql = db + "_" + save_date + ".sql"
        dumpsql = ("mysqldump --databases %s > %s" %(db, savesql))
        run = subprocess.Popen(dumpsql, stdout=subprocess.PIPE, shell=True)
        out, err = run.communicate()
        print(savesql)
# =============================================================================
# Run
# =============================================================================
if __name__ == '__main__':
#Create database, root must be set password 'read set_database_for_TAROT.txt'
    if False:
        import json
        Create_Database('Pydatabase')
        Create_Database('TCA')
        Create_Database('TCH')
        Create_Database('TRE')
    
#Create table    
        Create_Alert_User_Table()
        Create_Gws_Table()
        
        Create_Table(db='TCA', tb='gwevent')
        Create_Table(db='TCH', tb='gwevent')
        Create_Table(db='TRE', tb='gwevent')
        
    if False:        
#Call Database and Table
        Call_Database()
        Call_Table()
        Call_Table('TCA')
        
        Describe_Table(db='Pydatabase', tb='AlertUser')
        Describe_Table(db='Pydatabase', tb='GWs')
        Describe_Table(db='TCA', tb='gwevent')
        
#Insert Alert User
        Insert_Alert("This is the test", "S200105ae")

#Insert table Gws --> event.php        
        with open('/home/kanthanakorn/Documents/TSP/tarot/gws/gws_burst_time.txt', 'r') as rj:
            gws_dict = json.load(rj)
            
        for i,j in gws_dict.items():
            Insert_Gws([i, j['burst'].replace("T"," "), str(j['classification'][1:6]), 0, 0, 0])
            
#Insert table Ges --> event.php by event, go to /tarot/gws/and run
#        python3 gws_info.py S200105ae

    if False:
        pass
            
#    if False:
#        gws, gw = Call_GW_Event()
#        print(gws, gw)
#    

    if False:    
#        update alert user
        Update_Alert("This is the test", "S200105ae")

    if False:            
#        update tiles
#        update the rest of tiles , go to read note in fn.
#        [No.tca, No.tch, No.tre, superevent]
        S190412m = [149,0,100, 'S190412m']
        Update_Gws(S190412m)
        S190421ar = [0, 0, 72, 'S190421ar']
        Update_Gws(S190421ar)
        
    if False:            
#        update all tiles that TAROT has obs
        from dict_attributes import dict_tiles
        all_tiles = dict_tiles()
        for i in all_tiles:
            for m, n in i.items():
                print("Update tiles: %s" %m)
                Update_Gws(n)
                
    if False:
#        update tiles by import search_gw_alert_grandma.py
        from tarot.search.search_gw_alert_grandma import TAROT_Retrieve
        run = TAROT_Retrieve()
        run.Check_Log_link('TCA')
        run.Read_Log('TCA')
        
#        Read file fits produced by telescope
        tca_tile = []
        for i in run.table['sname']:
            tca_tile.append(i[0:9])
        
        run.Check_Log_link('TCH')
        run.Read_Log('TCH')
        tch_tile = []
        for i in run.table['sname']:
            tch_tile.append(i[0:9])
            
        run.Check_Log_link('TRE')
        run.Read_Log('TRE')
        tre_tile = []
        for i in run.table['sname']:
            tre_tile.append(i[0:9])
            
#        Read GW events from GWs.txt or GWs_spf.txt, choose what you need
        fileA = False
        fileB = False
        fileC = True
        fileD = False
        
        if fileA:
            with open("/home/kanthanakorn/Documents/TSP/Gws.txt") as rg:
                rgl = rg.read().splitlines()
        else:
            pass
    
        if fileB:
            with open("/home/tarot/TSP/Gws.txt") as rg:
                rgl = rg.read().splitlines()
        else:
            pass

#       Run some event
        if fileC:
            with open("/home/tarot/TSP/Gws_spf.txt") as rg:
                rgl = rg.read().splitlines()
        else:
            pass
        
        if fileD:
            with open("/home/kanthanakorn/Documents/TSP/Gws_spf.txt") as rg:
                rgl = rg.read().splitlines()
        else:
            pass
                
        for i in rgl:
            tca = tca_tile.count(i)
            tch = tch_tile.count(i)
            tre = tre_tile.count(i)
            tiles = [tca, tch, tre, str(i.replace("_", " "))]
            print("Update tiles")
            print(tiles)
            Update_Gws(tiles)
            
    if False:
#        update tiles in sys.argv[1]
        from tarot.search.search_gw_alert_grandma import TAROT_Retrieve
        import sys
        gw=sys.argv[1]
        # gw= 'S200115j'
        run_tca = TAROT_Retrieve()
        try:
            run_tca.Check_Log_link('TCA', gw)
            run_tca.Read_Log('TCA')
            tca_tile = len(run_tca.table['sname'])
        except:
            tca_tile = 0
            
        run_tch= TAROT_Retrieve()
        run_tch.Check_Log_link('TCH', gw)
        try:
            run_tch.Read_Log('TCH')
            tch_tile = len(run_tch.table['sname'])
        except:
            tch_tile = 0
            
        run_tre = TAROT_Retrieve()
        try:
            run_tre.Check_Log_link('TRE', gw)
            run_tre.Read_Log('TRE')
            tre_tile = len(run_tre.table['sname'])
        except:
            tre_tile = 0
        
        tiles = [tca_tile, tch_tile, tre_tile, gw]
        print("Update tiles")
        print(tiles)
        Update_Gws(tiles)
    
    if False:
#        update tiles from /tmp/*_FITS_files.txt and specific event 
        import sys
        from astropy.table import Table
        from tarot.search.search_gw_alert_grandma import TAROT_Retrieve
        try:
            tmp_tca = Table.read('/tmp/TCA_FITS_files.txt', format='ascii')
            tca_tile = len(tmp_tca)
        except(FileNotFoundError):
            tca_tile = 0
        
        try:
            tmp_tch = Table.read('/tmp/TCH_FITS_files.txt', format='ascii')
            tch_tile = len(tmp_tch)
        except(FileNotFoundError):
            tch_tile = 0
            
        try:
            tmp_tre = Table.read('/tmp/TRE_FITS_files.txt', format='ascii')
        except(FileNotFoundError):
            tre_tile = 0
            
        gw = 'S200115j'
        # gw = sys.argv[1]        
        
        tiles = [tca_tile, tch_tile, tre_tile, gw]
        print("Update tiles")
        print(tiles)
        # Update_Gws(tiles)
        
#    update gws that is not on the link
    if False:
        S190412m = [149,0,100, 'S190412m']
        Update_Gws(S190412m)
        S190421ar = [0, 0, 72, 'S190421ar']
        Update_Gws(S190421ar)

    if False:            
#        insert tiles
#        [superevent, date, [classification], No.tca, No.tch, No.tre]
#        file is in json 'gws_burst_time.txt' and 'gws_info.py in 'gws'
        S200101 = ['S200101', '2020-01-01 00:00:01', '[0, 0, 100, 0, 0]', 9, 9, 9]
        Insert_Gws(S200101)

    if False:
        Alter_Table_Name(tb="gwevent", new_tb="S200110a", db = "TCA")
        
        #then create 'gwevent for late data analysis
        Create_Table(db='TCA', tb='gwevent')
        
    if False:
        Alter_Table_Name(tb="gwevent", new_tb="S200110a", db = "TRE")
        
        #then create 'gwevent for late data analysis
        Create_Table(db='TRE', tb='gwevent')
