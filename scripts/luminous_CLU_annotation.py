"""
ZTF CLU fritz comment annotator that checks all the phootmetry from friz.scinece and posts a comment if the source is/or not a luminous CLU
event (given that ZTF CLU now is a volume-luminosity-limited SN Survey)
Author: Anastasios Tzanidakis (atzanida@caltech.edu) & all rights reserved California Institute of Technology
"""

import numpy as np
import os, glob
import requests, json
from lxml import html
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii

# Define global parameters
global TOKEN, BASEURL, USR_LAST_NAME

with open('user_info.json') as usr:
    # User Infomation
    usr_data = json.load(usr)

GETTOKEN = usr_data['FritzToken']
USR_LAST_NAME = usr_data['User_last_name']
BASEURL = 'https://fritz.science/'

def api(method, endpoint, data=None):
    ''' Info : Basic API query, takes input the method (eg. GET, POST, etc.), the endpoint (i.e. API url)
               and additional data for filtering
        Returns : response in json format
        CAUTION! : If the query doesn't go through, try putting the 'data' input in 'data' or 'params'
                    argument in requests.request call
    '''
    headers = {'Authorization': f'token {GETTOKEN}'}

    response = requests.request(method, endpoint, json=data, headers=headers)

    return response.json()

def get_source_api(ztfname, comments=False):
    ''' Info : Query a single source, takes input ZTF name
        comments (bool): If True it will also include the source comments
        Returns : all basic data of that source (excludes photometry and spectra,
                  includes redshift, classification, comments, etc.)
    '''
    if comments:
        url = BASEURL + 'api/sources/%s?includeComments=true'%ztfname
    else:
        url = BASEURL+'api/sources/%s'%ztfname

    response = api('GET', url)
    return response['data']

def get_source_coords(ztfname):
    """Generate ra and dec coordinates given ZTF source"""
    source = get_source_api(ztfname)
    return (source['ra'], source['dec'])

def post_source_comment(ztfname, comment):
    """ Info: Post comment to a specific ZTF name."""
    url = BASEURL + "api/comment"

    filt = {"obj_id": ztfname, 'text':comment}
    response = api('POST',url, filt)
    return (response)

def get_source_comment(number):
    """ Info: Post comment to a specific ZTF name."""
    url = BASEURL + "api/comment/" + number
    response = api('GET',url)
    return (response)

def get_source_photometry(ztfname, extrapolate_photometry=False, abs_mag_genator=False):
    """Fetch .json all photometry issued for a given ZTF source"""
    url = BASEURL + "api/sources/" + ztfname + "/photometry"
    response = api("GET", url)
    mjd, filt, mag, mag_err = [], [], [], []
    if extrapolate_photometry:
        phot = response['data']
        for alert in phot:
            if alert['mag']!=None:
                mjd.append(alert['mjd'])
                filt.append(alert['filter'])
                mag.append(alert['mag'])
                mag_err.append(alert['magerr'])

        # Convert to numpy arrays
        mjd, filt, mag, mag_err = np.array(mjd), np.array(filt), np.array(mag), np.array(mag_err)
        return (Table([mjd, filt, mag, mag_err], names=("mjd", "filter", 'mag', 'mag_err')))
    else:
        return response['data']

def fetch_redshift(ztfname):
    z = get_source_api(ztfname)['redshift']
    if not z:
        return (False)
    else:
        return (z)

def dist_mod_mag(app_mag, distance):
    """ Calculate the distance modulus of cosmological sources """
    return (app_mag - 5*np.log10(distance)-25)

def redshift_to_distance(z):
    # Define cosmological constanst from Astropy
    # H0 Estimates used on fritz.science
    cosmo = FlatLambdaCDM(H0=67.7 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    return cosmo.luminosity_distance(z).value

def CLU_luminous(ztfname, diagnostic_plot=True):
    """Return bool.

    If Mabs>=-17 (False)
    If Mabs<-17 (True)
    Show diagnotic: True (show the app mag and abs mag plot! with the selected point!)
    """
    # Fetch photometry & redshift
    photometry, source_z = get_source_photometry(ztfname, extrapolate_photometry=True), fetch_redshift(ztfname)

    if source_z!=False:
        # Convert magnitudes to absolute magnitude
        lum_dist = redshift_to_distance(source_z)
        Mabs = dist_mod_mag(photometry['mag'].data, lum_dist)

        # Lightcurve parameters
        unq_filters = np.unique(photometry['filter'].data)
        filt = photometry['filter'].data
        Mapp = photometry['mag'].data
        Mapp_err = (photometry['mag_err'].data)
        phase = photometry['mjd'].data - photometry['mjd'].data[0]

        # Select peak in lightcurve (first maximum)
        peak = np.argmin(Mabs)

        if diagnostic_plot:
            plt.figure(figsize=(6,5))
            for un in unq_filters:
                if un=='ztfr':
                    color='red'
                elif un=='ztfg':
                    color='green'
                else:
                    color='darkorange'
                plt.errorbar(phase[filt==un]-phase[peak], Mabs[filt==un], yerr=Mapp_err[filt==un], fmt='o',
                            label='%s'%un, color=color, mec='k', capsize=3)
            plt.legend(bbox_to_anchor=(1.4,1))
            plt.ylim(plt.ylim()[::-1])
            plt.ylabel(r"$M_{abs}$ (AB)", fontsize=16)
            plt.xlabel(r"Phase [days]", fontsize=16)
            plt.axvline(0, color='orange', lw=5, alpha=0.5)
            plt.xlim(-30, 30)
        else:
            if Mabs[peak]<=-17:
                return (True)
            else:
                return (False)
    else:
        return (False)

def get_all_sources(program_id):
    """ Info: Post comment to a specific ZTF name."""
    url = BASEURL + 'api/sources?group_ids=%s&saveSummary=true'%program_id
    print (url)
    response = api('GET',url)
    return (response)

def update_source_comment(comment_id, comment, obj_id, author_id):
    url = BASEURL + "api/comment/%s"%comment_id
    filt = {"text": comment, 'author_id':author_id, "obj_id":obj_id} # author_ID:37 (Andy Tzanidakis)
    response = api('PUT',url, filt)
    return (response)

def contain_luminosity_comment(ztfname):
    """
    Given a ZTF_Id if the "[Luminous_CLU_Event]" comment is identified it will return the comment_id. If not identified, it will return a bool (False)
    """
    source = get_source_api(ztfname, comments=True)
    comm_id = []
    for s in source['comments']:
        if s['text']=="[Luminous_CLU_Event] False" or s['text']=="[Luminous_CLU_Event] True":
            comm_id.append(s['id'])
    if len(comm_id)>0:
        return (comm_id[0])
    else:
        return (False)

def get_user_id(last_name, group_id=43):
    """
    Given the users last name and group_ID. This function will return the user_id specific to the last name.
    Params:
    last_name (str): Last name of the user (as it appears on fritz.science)
    group_id (float): Group ID (i.e 43 is Caltech Census of the Local Universe)
    """
    url = BASEURL + 'api/groups/%s'%group_id
    response = api("GET", url)
    usr_id = []
    for user in response['data']['users']:
        if user['last_name']==last_name:
            usr_id.append(user['id'])
    if len(usr_id)>0:
        return (usr_id)
    else:
        return (False)

def main():

    usr_id = get_user_id(USR_LAST_NAME, group_id=43)
    # Fetch all CLU sources
    all_clu_sources = get_all_sources(43)

    # Define Start Date you want to query
    t_start_query  = Time("2021-01-01")

    for source in all_clu_sources['data']['sources']:
        T = Time(source['saved_at']) # as astropy.time

        if T>=t_start_query: # if this is greater or equal to the date
            cand_id, date_id = source['obj_id'], source['saved_at']

            # Check photometry and return boolean if luminous or not (i.e TRUE/FALSE)
            lum_clu_check = CLU_luminous(cand_id, diagnostic_plot=False)

            # Check if Luminosity CLU comment exist
            c_id = contain_luminosity_comment(cand_id)

            if c_id!=False: # found a comment
                # Update this comment with it's given comment_id
                update_comment = update_source_comment(c_id, "[Luminous_CLU_Event] %s"%lum_clu_check, cand_id, usr_id)
                print ("Updated comment for %s"%cand_id)
            else:
                # Post new comment
                comment_poster = post_source_comment(cand_id, "[Luminous_CLU_Event] %s"%lum_clu_check)

            # Summary of data posting!
            print ("#########")
            print ("Checking luminosity of source ", cand_id, "stored at ", date_id)
            print ("Luminosity Check:", lum_clu_check)

if __name__ == "__main__":
    main()
