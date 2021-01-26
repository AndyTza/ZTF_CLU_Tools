import numpy as np
import os, glob
import requests, json
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy import constants as const
from astropy.table import Table
import threading
import time
from tqdm import tqdm
import warnings
from datetime import date

today = date.today()
warnings.filterwarnings('ignore')

# Define global parameters
global TOKEN, BASEURL, USR_LAST_NAME

with open('user_info.json') as usr:
    # User Infomation
    usr_data = json.load(usr)

GETTOKEN = usr_data['user']['FritzToken']
USR_LAST_NAME = usr_data['user']['user_last_name']
BASEURL = 'https://fritz.science/'

global obj_id, saved_date, ra, dec, z, clu_d, classification, luminous_event, peak_app_mag, classification_prob, peak_abs_mag

# List of Empty Parameters
obj_id, saved_date, ra, dec = [], [], [], []
z, clu_d = [], []
classification, classification_prob = [], []
luminous_event, peak_app_mag, peak_abs_mag = [], [], []

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

def get_all_sources(program_id):
    """
    Fetch all the obj_id and save dates of all sources for a given program_id.
    program_id (int): program_id of your program (i.e 43 for Caltech Census of the Local Univese)
    """
    url = BASEURL + f'api/sources?group_ids={program_id}&saveSummary=true'
    print (f"Query API URL: {url}")
    response = api('GET',url)
    return (response)

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

def get_user_id(last_name, group_id=43):
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

def CLU_luminous(ztfname, source_z, luminosity_threshold=-17.0):
    """
    This function will return a bool (True/False) according a to a luminosity threshold.
    The estimates of the luminosity take into account the Peak absolute magnitude (peak light) + error.

    Input
    -----
    ztfname (str): obj_id of the ZTF source
    soure_z (float): redshift of the ZTF source
    luminosity_threshold (float): luminosity threshold you want to test your peak absolute magnitude

    Output:
    luminous event (bool): True if it's more luminous than the luminosity threshold & False if fainter of the luminosity threshold
    """
    # Fech photometry & redshift
    photometry = get_source_photometry(ztfname, extrapolate_photometry=True)

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
        phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

        # Select peak in lightcurve (first maximum)
        peak = np.argmin(Mabs)
        peak_app = np.argmin(Mapp) # peak in apparent magnitude space

        # Check if source is luminous by considering the absolute magnitude + error
        is_luminous_error = (Mabs[peak] + Mapp_err[peak])<=luminosity_threshold

        return (is_luminous_error, Mabs[peak])
    else:
        return (False, None)

def peak_mag_CLU(ztfname):
    """
    Return the peak apparent magnitude given the ztf name. This will return the first maximum detection.
    Input
    -----
    ztfname (str): Obj_id of a given source

    Output
    ------
    peak_app_mag (float): peak apparent magnitude
    """
    # Query photometry
    photometry = get_source_photometry(ztfname, extrapolate_photometry=True)

    try:
        nan_index = photometry['mag'].data!=None # select index where no None's

        # Lightcurve parameters
        unq_filters = np.unique(photometry['filter'].data[nan_index])
        filt = photometry['filter'].data[nan_index]
        Mapp = photometry['mag'].data[nan_index]
        Mapp_err = photometry['mag_err'].data[nan_index]
        phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

        peak_app = np.argmin(Mapp) # peak in apparent magnitude space

        return Mapp[peak_app]
    except:
        return (None)

def query_CLU(ZTF_name, saved_at):
    """
    This function will append the parameters of a given ZTF souce to the
    global parameters (see query output).

    Query Input:
    -----------
    ZTF_name (str): obj_id of the given ZTF source
    saved_at (str): saved date of the given source

    Query Output:
    -------------
    obj_id (str): obj_id of a given ZTF source
    saved_at (str): date saved on fritz.scince
    ra, dec (float): RA, DEC coordinates of the source (mean position)
    z (float): redshift of CLU host galaxy
    clu_d (float): sky seperation of the transient with respect to the CLU host (arcsec)
    classification (str): the most recent classification on fritz.science
    classification_prob(float): probability of the most recent classification
    peak_app_mag (float): peak apparent magnitude of the lightcurve
    peak_abs_mag (float): peak absolute magnitude of the ligtcurve given the most recent redshift estimate
    luminous_event (str): if source is luminous or subluminous given a -17 luminosity cut
    """
    source = get_source_api(ZTF_name, comments=False) # Fetch source infomation

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

    if len(clu_z)>0:
        # Make into numpy arrays
        clu_z, clu_d_gal_sep = np.array(clu_z), np.array(clu_d_gal_sep)

        # Find the indicies of the minimum distance(seperation in arcsec)
        find_nearest = np.argmin(clu_d_gal_sep)

        z.append(clu_z[find_nearest])
        clu_d.append(clu_d_gal_sep[find_nearest])

        # Check if luminous
        luminosity_info = CLU_luminous(ZTF_name, clu_z[find_nearest], luminosity_threshold=-17.0)
        luminous_event.append(luminosity_info[0])
        peak_abs_mag.append(luminosity_info[1])

    else:
        z.append(source['redshift'])
        clu_d.append(None)

        luminosity_info = CLU_luminous(ZTF_name, source['redshift'], luminosity_threshold=-17.0)

        # Cannot run luminosity check....
        luminous_event.append(luminosity_info[0])
        peak_abs_mag.append(luminosity_info[1])


    # Fetch the most *recent* classification
    clf = source['classifications']

    if clf:
        classification.append(clf[-1]['classification']) # classification
        classification_prob.append(clf[-1]['probability']) # classification probability 0-1
    else:
        classification.append(None)
        classification_prob.append(None)

    # Fetch peak apparent magnitudes
    peak_apparent_magnitude = peak_mag_CLU(ZTF_name)

    if peak_apparent_magnitude!=None:
        peak_app_mag.append(peak_apparent_magnitude)
    else:
        peak_app_mag.append(None)

    # Append Parameters
    obj_id.append(ZTF_name) # Obj_id of the ZTF source
    saved_date.append(saved_at) # Date it was saved from
    ra.append(source['ra']) # R.A coord
    dec.append(source['dec']) # DEC coord

def main():
    sources = get_all_sources(43) # fetch all sources for ZTF CLU experiment

    # All souces (names & dates) of all saved ZTF CLU candidates
    names = [s['obj_id'] for s in sources['data']['sources']]
    dates = [s['saved_at'] for s in sources['data']['sources']]

    list_threads = [] # list of threads for pooling
    buff = 0 # buffer time (NOTE: API sometimes crashes if request is to frequent)
    for i in tqdm(range(len(names))):
        buff +=1
        if buff>=50: # feed to pool 50 sources a time
            time.sleep(12) # sleep for ~12 seconds (for ~3000 sources this takes roughly 10 minutes)
            buff = 0 # restart buffer

        t = threading.Thread(target=query_CLU, args=[names[i], dates[i]])
        t.start()
        list_threads.append(t)

    for thread in list_threads:
        thread.join()

    # Extrapolate and store all data
    param_names = ("ZTF_id", "saved_date", "ra", "dec",
     "z", "clu_d_host", "classification", "classification_prob",
      "peak_app_mag", "peak_abs_mag", "luminous_event")

    data_table = Table([obj_id, saved_date, ra, dec, z, clu_d,
     classification,classification_prob, peak_app_mag, peak_abs_mag,
     luminous_event], names=param_names)

    data_table.write(f"CLU_QUERY_{today.year}_{today.month}_{today.day}.ascii", format='ascii', overwrite=True)

if __name__ == "__main__":
    main()
