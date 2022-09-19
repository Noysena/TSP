#TSP = Transient Search Package
------------------------------------------------------------------------------------------------
#Need python3-pip
pip3
pip3 install --user numpy 
pip3 install --user scipy 
pip3 install --user beautifulsoup4
pip3 install --user h5py 
pip3 install --user pyyaml 
pip3 install --user matplotlib 
pip3 install --user pytz 
pip3 install --user scikit-image 
pip3 install --user pandas
pip3 install --user objgraph 
pip3 install --user setuptools 
pip3 install --user mock 
pip3 install --user astropy
pip3 install --user lalsuite
pip3 install --user astroquery
pip3 install --user numba
pip3 install --user healpy
pip3 install --user cryptography 
pip3 install --user scikit-learn
pip3 install --user pyds9
pip3 install --user mysql-connector
pip3 install --user ligo.skymap
------------------------------------------------------------------------------------------------
#Before install ligo.skymap
sudo apt-get install libgsl23 libgsl-dbg libgsl-dev
sudo apt-get install astrometry.net

pip3 show ligo.skymap
#Then look for ligo/skymap
cp ~tarot/gws/ligo_area.py /home/kanthanakorn/.local/lib/python3.7/site-packages/ligo/skymap
cp ~tarot/gws/ligo_distance.py /home/kanthanakorn/.local/lib/python3.7/site-packages/ligo/skymap
------------------------------------------------------------------------------------------------
#To test package
$python3 setup.py test

------------------------------------------------------------------------------------------------
#To install package
$pip3 install --user -e .

------------------------------------------------------------------------------------------------
#Install XAMPP --> https://www.apachefriends.org/index.html
#file location --> 
sudo apt install net-tools
sudo apt-get install mariadb-client-core-10.3
#set password for sql --> https://www.youtube.com/watch?v=ZSMRmvIUhdE
#set password 'root' 3 of them

#go to --> /opt/lampp/phpmyadmin/config.inc.php
#edit
# --> $cfg['Servers'][$i]['password'] = '';
# --> $cfg['Servers'][$i]['password'] = 'deckanpolar2!';

#start --> sudo /opt/lampp/lampp start
#stop --> sudo /opt/lampp/lampp stop

#test XAMPP --> http://localhost
#put folder of web here
/opt/lampp/htdocs

#Connect to sql
mysql -h 127.0.0.1 -P 3306 -u root -p

#set permissions to htdocs
drwxr-xr-x  5 root   root      4096 Jan 29 17:20 htdocs
#To
drwxr-xr-x  5 kanthanakorn kanthanakorn    4096 Jan 29 17:20 htdocs

sudo chown -R kanthanakorn:kanthanakorn /opt/lammp/htdocs
#edit /opt/lampp/etc/httpd.conf
sudo gedit /opt/lampp/etc/httpd.conf
#find and replace with name you change: line 173-174
User daemon
Group daemon
#To
User kanthanakorn
Group kanthanakorn

######### Run tarot TSP_SQL.py #############
#Create database and table for website.
Create_Database('Pydatabase')
Create_Database('TCA')
Create_Database('TCH')
Create_Database('TRE')
Create_Table(db='TCA', tb='gwevent')
Create_Table(db='TCH', tb='gwevent')
Create_Table(db='TRE', tb='gwevent')

#Create table GWs table for event.php
Create_Gws_Table(tb="GWs")
Create_Alert_User_Table(tb="AlertUser")
#--> Insert_Alert("Alert", "S202020") , but Update_Alert("", "") when you have had it.
------------------------------------------------------------------------------------------------

#To run TSP, simply do...
import tarot
print(tarot.TAROT_PIP.MAIN())

------------------------------------------------------------------------------------------------

#incase no class just simple do...
print(tarot.MAIN())
#PATH is not success to run directly form terminal yet!

