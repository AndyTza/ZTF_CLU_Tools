from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import bs4 as bs
import re
from selenium import webdriver
import time
from astropy.io import ascii
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm
from astropy.table import Table


CLU = # Load CLU list of candidates

def redshift_to_distance(z):
    # Define cosmological constanst from Astropy
    cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    return (cosmo.luminosity_distance(z))


def HST_img_query(ra, dec, radius):
    """
    Return the url querying the HST Arhcive Image website.

    Input:
    -------
    RA (J2000), DEC (J2000)
    radius (arcsec): float

    Output:
    ------
    bool: if True, it found result in images on HLA & will return the url to the images
          if False, it did not find images on HLA
    """

    # Convert radius into degrees
    base_url = "https://hla.stsci.edu/hlaview.html#Images|filterText%3D%24filterTypes%3D|posfilename=&poslocalname=&posfilecount=&listdelimiter=whitespace&listformat=degrees&"
    coord_url = f"&RA={ra}&Dec={dec}&Radius={radius*u.arcsec.to(u.deg).value}&"
    inst_url = "inst-control=all&inst=ACS&inst=ACSGrism&inst=WFC3&inst=WFPC2&inst=NICMOS&inst=NICGRISM&inst=COS&inst=WFPC2-PC&inst=STIS&inst=FOS&inst=GHRS&imagetype=best&prop_id=&"
    data_url = "spectral_elt=&proprietary=both&preview=1&output_size=256&cutout_size=12.8|ra=&dec=&sr=&level=&image=&inst=ACS%2CACSGrism%2CWFC3%2CWFPC2%2CNICMOS%2CNICGRISM%2CCOS%2CWFPC2-PC%2CSTIS%2CFOS%2CGHRS&ds="

    url_new = base_url + coord_url + inst_url + data_url

    driver = webdriver.Chrome(executable_path='') # add your crhomedriver path <HERE>
    driver.get(url_new)
    time.sleep(2) # snooz

    html = driver.page_source
    soup = bs.BeautifulSoup(html)

    k = soup.find("div", {'id': 'output'})

    found = False
    ent = []
    for j in k.text.split("\t"):
        op = j.split(" ")
        if len(op)==2:
            if op[0]=='Results':
                ent.append(op[1])

    ent = ent[0]
    if float(ent[0])==0:
        found = False
        #print (float(ent[0]))
    else:
        found = True

    driver.quit()
    return (found, url_new)


p0_name, p0_clf = [], []
p0_url = []

for i in tqdm(range(len(CLU_near))):
    s = CLU_near[i]
    hst = HST_img_query(s['ra'], s['dec'], 30)

    if hst[0]:
        p0_name.append(s['ZTF_id'])
        p0_clf.append(s['classification'])
        p0_url.append(hst[1])

T = Table([p0_name, p0_clf, p0_url])
T.write('table_HST.ascii', format='ascii')
