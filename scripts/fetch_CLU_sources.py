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
global TOKEN, BASEURL, USR_LAST_NAME, db

with open('user_info.json') as usr:
    # User Infomation
    usr_data = json.load(usr)

GETTOKEN = usr_data['user']['FritzToken']
USR_LAST_NAME = usr_data['user']['user_last_name']
BASEURL = 'https://fritz.science/'

# Empty dict.
db = {}
db['source'] = []

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

        if len(mag):
            return (Table([mjd, filt, mag, mag_err], names=("mjd", "filter", 'mag', 'mag_err')))
        else:
            return (False)
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
        #phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

        if len(filt):
            # Select peak in lightcurve (first maximum)
            peak = np.argmin(Mabs)

            # Check if source is luminous by considering the absolute magnitude + error
            is_luminous_error = (Mabs[peak] + Mapp_err[peak])<=luminosity_threshold

            return (is_luminous_error, Mabs[peak])
        else:
            return (False, None)
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
    if photometry: # if photometry table is not empty
        try:
            nan_index = photometry['mag'].data!=None # select index where no None's

            # Lightcurve parameters
            unq_filters = np.unique(photometry['filter'].data[nan_index])
            filt = photometry['filter'].data[nan_index]
            Mapp = photometry['mag'].data[nan_index]
            Mapp_err = photometry['mag_err'].data[nan_index]
            phase = photometry['mjd'].data[nan_index] - photometry['mjd'].data[nan_index][0]

            # Try to get the peak magnitude in r-band, if r-band does not exist, any other photometry is okay
            contain_red, contain_green, contain_i = False, False, False
            for f in filt:
                if f=="ztfr": # rband
                    contain_red = True
                elif f=='ztfg':
                    contain_green = True
                elif f=='ztfi':
                    contain_i = True

            if contain_red==True:
                Mapp_phot_band = Mapp[filt=='ztfr']
            elif contain_red==False and contain_green==True:
                Mapp_phot_band = Mapp[filt=='ztfg']
            elif contain_red==False and contain_green==False and contain_i==True:
                Mapp_phot_band = Mapp[filt=='ztfi']

            peak_app = np.argmin(Mapp_phot_band) # peak in apparent magnitude space

            return Mapp_phot_band[peak_app]
        except:
            return (None)
    else:
        return (False)

def query_CLU_dict(ZTF_name, saved_at):
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

        # Check if luminous
        luminosity_info = CLU_luminous(ZTF_name, clu_z[find_nearest], luminosity_threshold=-17.0)

        # Assign parameters
        param_z = float(clu_z[find_nearest])
        param_clu_d = float(clu_d_gal_sep[find_nearest])
        param_luminous_event = luminosity_info[0]
        param_peak_abs_mag = luminosity_info[1]

    else:
        luminosity_info = CLU_luminous(ZTF_name, source['redshift'], luminosity_threshold=-17.0)

        # Assign parameters
        param_z = float(source['redshift'])
        param_clu_d = None
        param_luminous_event = luminosity_info[0]
        param_peak_abs_mag = luminosity_info[1]

    # Fetch peak apparent magnitudes
    peak_apparent_magnitude = peak_mag_CLU(ZTF_name)

    if peak_apparent_magnitude!=None:
        param_peak_app_mag = peak_apparent_magnitude
    else:
        param_peak_app_mag = None

    # Fetch the most *recent* classification
    clf = source['classifications']

    if clf:
        param_classification = clf[-1]['classification']
        param_prob = clf[-1]['probability']
    else:
        param_classification = None
        param_prob = None

    # Assign last parameters
    param_obj_id = ZTF_name
    param_saved_date = saved_at
    param_ra = source['ra']
    param_dec = source['dec']

    # Extrapolate to main dict.
    db['source'].append({"ZTF": param_obj_id, "saved_date": param_saved_date,
                         "ra":param_ra, "dec":param_dec,
                        "z": param_z, "clu_d": param_clu_d,
                         "classification": param_classification, "probability": param_prob,
                        "peak_abs_mag": param_peak_abs_mag,
                         "luminous_event": str(param_luminous_event), "peak_app_mag": param_peak_app_mag})

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
            time.sleep(13) # sleep for ~10 seconds (for ~3000 sources this takes roughly 10 minutes)
            buff = 0 # restart buffer

        t = threading.Thread(target=query_CLU_dict, args=[names[i], dates[i]])
        t.start()
        list_threads.append(t)

    for thread in list_threads:
        thread.join()

    # Extrapolate and store all data
    param_names = ("ZTF_id", "saved_date", "ra", "dec",
     "z", "clu_d_host", "classification", "classification_prob",
      "peak_app_mag", "peak_abs_mag", "luminous_event")

    # Extrapolate from the dict. to list
    obj_id = [s['ZTF'] for s in db['source']]
    saved_date = [s['saved_date'] for s in db['source']]
    ra = [s['ra'] for s in db['source']]
    dec = [s['dec'] for s in db['source']]
    z = [s['z'] for s in db['source']]
    clu_d = [s['clu_d'] for s in db['source']]
    classification = [s['classification'] for s in db['source']]
    classification_prob = [s['probability'] for s in db['source']]
    peak_app_mag = [s['peak_app_mag'] for s in db['source']]
    peak_abs_mag = [s['peak_abs_mag'] for s in db['source']]
    luminous_event = [s['luminous_event'] for s in db['source']]

    data_table = Table([obj_id, saved_date, ra, dec, z, clu_d,
     classification, classification_prob, peak_app_mag, peak_abs_mag,
     luminous_event], names=param_names)

    data_table.write(f"CLU_QUERY_{today.year}_{today.month}_{today.day}.ascii", format='ascii', overwrite=True)

if __name__ == "__main__":
    main()
