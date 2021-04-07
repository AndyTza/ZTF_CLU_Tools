#!/usr/bin/env python
"""
Query the most recent ZTF CLU saved sources & check if they pass a luminosity threshold the default is -17.
This code will leave an annotation in the format: [Luminous_CLU_Event]: True/False

Author: Anastasios Tzanidakis (atzanida@caltech.edu)
"""
import numpy as np
import os, glob
import requests, json
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from astropy.table import Table
from astropy.time import Time
from astropy.io import ascii
import argparse

# Define global parameters
global TOKEN, BASEURL, USR_LAST_NAME

with open('user_info.json') as usr:
    # User Infomation
    usr_data = json.load(usr)

GETTOKEN = usr_data['user']['FritzToken']
USR_LAST_NAME = usr_data['user']['user_last_name']
BASEURL = 'https://fritz.science/'

parser = argparse.ArgumentParser(description="Luminous ZTF CLU Annotator")
parser.add_argument('-d', '--date', type=str, help='Add the date of your query [YYYY-MM-DD]')
args = parser.parse_args()

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
        Returns : all basic data of that source
    '''
    if comments:
        url = BASEURL + f'api/sources/{ztfname}?includeComments=true'
    else:
        url = BASEURL+f'api/sources/{ztfname}'

    response = api('GET', url)
    return response['data']

def post_source_comment(ztfname, comment):
    """ Info: Post comment to a specific ZTF name."""
    url = BASEURL + "api/comment"

    filt = {"obj_id": ztfname, 'text':comment}
    response = api('POST',url, filt)
    return (response)

def get_source_photometry(ztfname, extrapolate_photometry=False):
    """ Fetch the photometry issued by ZTF for a given ZTF source.
    Params:
    -------
    ztfname (str): obj_id of the given ZTF source
    extrapolate_photometry (bool): Extrapolate the photometry from the .json file and convert it to astropy.Table format

    Output:
    ------
    astropy.Table containing mjd, filter, mag, mag_err (apparent magnitude), if extrapolate_photometry: True
    .json containing all magnitude information (including upper limits), if extrapolate_photometry: False
    """
    url = BASEURL + "api/sources/" + ztfname + "/photometry" # base-url to fetch data
    response = api("GET", url)

    if extrapolate_photometry:
        phot = response['data'] # .json file with all the photometric information
        mjd = [alert['mjd'] for alert in phot if alert['mag']!=False]
        filt = [alert['filter'] for alert in phot if alert['mag']!=False]
        mag = [alert['mag'] for alert in phot if alert['mag']!=False]
        mag_err = [alert['magerr'] for alert in phot if alert['mag']!=False]

        return (Table([mjd, filt, mag, mag_err], names=("mjd", "filter", 'mag', 'mag_err')))
    else:
        return (response['data'])

def fetch_redshift(ztfname):
    """
    Fetch the redshift of the nearest CLU galaxy (via CLU-auto-annoation). If not, will query the current redshift entry on Fritz.science

    Input:
    -----
    ztfname (str): obj_id of given ZTF source

    Output
    -----
    redshift (float): If redshift has been found on on fritz.science - or False (bool) if no redshift is found
    """

    source = get_source_api(ztfname, comments=False) # Fetch source infomation

    # Extrapolate CLU Redshift with the nearest crossmatch
    source_annotation = source['annotations'] # Fetch annotations
    clu_z, clu_d_gal_sep = [], []
    for annotation in source_annotation:
        annotation_data = annotation['data']
        for key in annotation_data.keys():
            if key=="cross_matches_CLU":
                for CLUdata in annotation_data['cross_matches_CLU']:
                    try:
                        clu_z.append(CLUdata['z'])
                        clu_d_gal_sep.append((CLUdata['coordinates'])['distance_arcsec'])
                    except:
                        if source['redshift']:
                            clu_z.append(source['redshift'])
                            clu_d_gal_sep.append(False)
                        else:
                            clu_z.append(False)
                            clu_d_gal_sep.append(False)

    # If the list has found redshift and seperation
    if len(clu_z)>0:
        clu_z, clu_d_gal_sep = np.array(clu_z), np.array(clu_d_gal_sep)

        # Nearest cross-match to CLU
        find_nearest = np.argmin(clu_d_gal_sep)

        return (clu_z[find_nearest])

    else:
        z = get_source_api(ztfname)['redshift'] # fetch source redshift
        if z:
            return (z)
        else:
            return (False)

def dist_mod_mag(app_mag, distance):
    """
    Calculate the absolute magnitude via the distance modulus.

    Input:
    -----
    app_mag (float): apparent magnitude of the source
    distance (float): distance of the source in Mpc

    Output:
    ------
    abs_mag (float): absolute magnitude
    """
    return (app_mag - 5*np.log10(distance)-25)

def redshift_to_distance(z):
    """
    Convert a given redshift to distance in Mpc, assuming H0=67.7 km/s/Mpc, and T_cmb=2.725 K
    Current cosmological constants used are similar to those on fritz.science

    Input:
    -----
    z (float): redshift of the source

    Output:
    ------
    distance (float): distance in Mpc
    """
    cosmo = FlatLambdaCDM(H0=67.7 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
    return cosmo.luminosity_distance(z).value

def CLU_luminous(ztfname, luminosity_threshold=-17.0):
    """
    This function will return a bool (True/False) according a to a luminosity threshold.
    The estimates of the luminosity take into account the Peak absolute magnitude (peak light) + error.

    Input
    -----
    ztfname (str): obj_id of the ZTF source
    luminosity_threshold (float): luminosity threshold you want to test your peak absolute magnitude

    Output:
    luminous event (bool): True if it's more luminous than the luminosity threshold & False if fainter of the luminosity threshold
    """
    # Fech photometry & redshift
    photometry = get_source_photometry(ztfname, extrapolate_photometry=True)
    source_z = fetch_redshift(ztfname)

    if source_z:
        # Convert magnitudes to absolute magnitude
        lum_dist = redshift_to_distance(source_z)

        nan_index = photometry['mag'].data!=None # select index where no None's

        # Fetch absolute magnitude of the source
        Mabs = dist_mod_mag(photometry['mag'].data[nan_index], lum_dist)

        # Lightcurve parameters
        unq_filters = np.unique(photometry['filter'].data[nan_index])
        filt = photometry['filter'].data[nan_index]
        Mapp = photometry['mag'].data[nan_index]
        Mapp_err = photometry['mag_err'].data[nan_index]

        try:
            if len(photometry['mjd'].data[nan_index])>0:
                phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]
            else:
                phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index] # likely only one detection
        except:
            return (False) # it failed at this step...

        try:
            # Select peak in lightcurve (first maximum)
            peak = np.argmin(Mabs)
        except:
            return (False)

        # Check if source is luminous by considering the absolute magnitude + error
        is_luminous_error = (Mabs[peak] + Mapp_err[peak])<=luminosity_threshold

        return (is_luminous_error)
    else:
        return (False)

def get_all_sources(program_id):
    """
    Fetch all the obj_id and save dates of all sources for a given program_id.
    program_id (int): program_id of your program (i.e 43 for Caltech Census of the Local Univese)
    """
    url = BASEURL + f'api/sources?group_ids={program_id}&saveSummary=true'
    print (f"Query API URL: {url}")
    response = api('GET',url)
    return (response)

def update_source_comment(comment_id, comment, obj_id, author_id):
    url = BASEURL + f"api/comment/{comment_id}"
    filt = {"text": comment, 'author_id':author_id, "obj_id":obj_id} # author_ID:37 (Andy Tzanidakis)
    response = api('PUT',url, filt)
    return (response)

def contain_luminosity_comment(ztfname):
    """
    Given a obj_id if the "[Luminous_CLU_Event]" comment is identified it will return the comment_id.
    If not identified, it will return a bool (False)

    Input:
    -----
    ztfname (str): obj_id of given ZTF source

    Output:
    ------
    contain_comment (bool): returns the bool if the [Luminous_CLU_Event] comment is found
    """
    source = get_source_api(ztfname, comments=True) # Fetch source with all the comments included

    # Check and append the comment id if found!
    comm_id = [s['id'] for s in source['comments'] if (s['text']=="[Luminous_CLU_Event] False") or (s['text']=="[Luminous_CLU_Event] True")]

    if comm_id:
        return (comm_id[0])
    else:
        return (False)

def get_user_id(last_name, group_id=43):
    """
    Given the users last name and group_id.
    This function will return the user_id specific to the last name.
    If last name is not found in the group, False (bool) will be returned

    Input:
    -----
    last_name (str): Last name of the user (as it appears on fritz.science)
    group_id (float): Group ID (i.e 43 is Caltech Census of the Local Universe)

    Output:
    -----
    usr_id (int): user_id unique to that user or False(bool) if not found
    """
    url = BASEURL + f'api/groups/{group_id}'
    response = api("GET", url)
    usr_id = [user['id'] for user in response['data']['users'] if user['last_name']==last_name]

    if usr_id:
        return (usr_id)
    else:
        return (False)

def main(date_request):
    # Fetch the user ID
    usr_id = get_user_id(USR_LAST_NAME, group_id=43)
    print (f"This is my user id: {usr_id}")

    # Fetch all CLU sources
    all_clu_sources = get_all_sources(43)

    # Define Start Date you want to query
    t_start_query  = Time(date_request, format='iso')
    counter = 0 # Counter for a buffer time
    for source in tqdm(all_clu_sources['data']['sources']):
        T = Time(source['saved_at']) # as astropy.time

        if T>=t_start_query: # if this is greater or equal to the date
            cand_id, date_id = source['obj_id'], source['saved_at']
            counter += 1

            if counter>=30:
                counter = 0 # restart counter...
                time.sleep(10) # sleep for 10 seconds...
            print (f"Checking: {cand_id}")
            # Check photometry and return boolean if luminous or not (i.e TRUE/FALSE)
            lum_clu_check = CLU_luminous(cand_id, luminosity_threshold=-17.0)

            # Check if Luminosity CLU comment exist
            c_id = contain_luminosity_comment(cand_id)

            if c_id!=False: # found a comment
                # Update this comment with it's given comment_id
                update_comment = update_source_comment(c_id, f"[Luminous_CLU_Event] {lum_clu_check}", cand_id, usr_id[0])
                print (f"comment status was: {update_comment['status']}")
                print (f"Updated comment for {cand_id} to {lum_clu_check}")
                print (" ")

            else:
                # Post new comment
                comment_poster = post_source_comment(cand_id, f"[Luminous_CLU_Event] {lum_clu_check}")

            # Summary of data posting!
            print ("Checking luminosity of source ", cand_id, "stored at ", date_id)
            print ("Luminosity Check:", lum_clu_check)
            print (f"fritz.science URL: https://fritz.science/source/{cand_id}")
            print (" ")

if __name__ == "__main__":
    main(args.date)
