#!/usr/bin/env python
"""
Upload a .txt file given a list of ZTF_ids to fritz.

Author: Anastasios Tzanidakis, Igor Andreoni, Yashvi Sharma
California Institute of Technology
"""
import numpy as np
import os, glob
import requests, json
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

parser = argparse.ArgumentParser(description="Upload Sources to Frtiz")
parser.add_argument('-g', '--group', type=int, help='add the group number you want to upload this targer (i.e 43)')
parser.add_argument('-f', '--file', type=str, help='name of the file that contains the ZTF objexts (note: make sure it is .txt/.ascii)')
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

def upload_source(ztfname, group_ids):
    ''' Info : Save an alert to given group(s), takes the ZTF name and list of group names
        Returns : response with status
    '''
    url = BASEURL+'api/alerts/ztf/'+ztfname
    filt = {'group_ids':[group_ids]}
    headers = {'Authorization': f'token {GETTOKEN}'}
    response = requests.request('POST', url, json=filt, headers=headers)
    return response.json()

def main(group_id, file_name):
    for source in ascii.read(file_name, format='fixed_width_no_header')['col1']:
        # First check if source is already saved to fritz
        if get_source_api(source):
            print (f"Sorry, it looks like {source} already exists on Frtiz.")
            continue
        else:
            upload_source(source, group_id)
            print (f"Congrats, source {source} has been uploaded to friz. Here is the url: https://fritz.science/source/{source}")

if __name__ == "__main__":
    main(args.group, args.file)
