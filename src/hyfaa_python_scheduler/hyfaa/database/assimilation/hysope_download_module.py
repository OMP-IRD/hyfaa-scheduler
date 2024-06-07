#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
#  This code has been developped by Magellium SAS
#
#  Licensing:
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>
#######################################################################

import requests, json
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import urllib
import io

from hyfaa.utils.fileenv import fileenv

hydroweb_credentials = {
    'user': fileenv('HYDROWEB_USER'),
    'pwd' : fileenv('HYDROWEB_PASSWORD'),
}
hydroweb_baserequest = 'http://hydroweb.theia-land.fr/hydroweb/authdownload?products='

hydroweb_next_service_url = "https://hydroweb.next.theia-land.fr/geoserver/HYDROWEB/ows?"

def get_hysope_wsh(sv_names, min_date=None, max_date=None, verbose=1):

    print('using hydroweb username {}'.format(hydroweb_credentials['user']))
    dico_out = dict()
    tstart0 = datetime.utcnow()
    for sv_name in set(sv_names):
        tstart = datetime.utcnow()
        req_txt = hydroweb_baserequest + sv_name + '&format=json'
        if min_date is not None:
            req_txt += '&sdate=%s'%min_date.strftime('%Y-%m-%d')
        if max_date is not None:
            req_txt += '&edate=%s'%max_date.strftime('%Y-%m-%d')
        req_txt += '&user=%s&pwd=%s'%(hydroweb_credentials['user'], hydroweb_credentials['pwd'])
        if verbose >= 1:
            print('URL {}'.format(req_txt))
        data_loc = requests.get(req_txt).json()['data']
        datetimes_loc = [datetime.strptime(el['date'] + 'T' + el['time'], '%Y-%m-%dT%H:%M') for el in data_loc]
        values_loc = [float(el['orthometric_height_of_water_surface_at_reference_position']) if el['orthometric_height_of_water_surface_at_reference_position'] is not None else None  for el in data_loc]
        uncertainties_loc = [float(el['associated_uncertainty']) if el['associated_uncertainty'] is not None else None for el in data_loc]

        #select only non NAN values
        ids_select = [ii for ii in range(len(datetimes_loc)) if (values_loc[ii] is not None) or (uncertainties_loc[ii] is not None)]
        datetimes_loc = [datetimes_loc[ii] for ii in ids_select]
        values_loc = [values_loc[ii] for ii in ids_select]
        uncertainties_loc = [uncertainties_loc[ii] for ii in ids_select]

        dico_out[sv_name] = {'dates': datetimes_loc, 'wsh': values_loc, 'wsh_uncertainty': uncertainties_loc}
        if verbose >= 1:
            print('Downloaded assimilation data from SV %s in %s'%(sv_name, datetime.utcnow()-tstart))
    if verbose >= 1:
        print('Downloaded all assimilation data in %s'%(datetime.utcnow()-tstart0))
    return dico_out


def get_hysope_wsh_hydrowebnext(sv_names, min_date=None, max_date=None, verbose=1):
    """
    Get water height from Hydroweb.next services. Should replace get_hysope_wsh function, that was using the old
    hydroweb API.

    The official API does not really provide the expected extraction
    services, so this is using the underlying geoserver WFS services, with ad-hoc CQL filter
    This will allow us to retrieve the full set of data in one single query

    :param sv_names: list of virtual stations identifiers
    :param min_date: datetime object
    :param max_date: datetime object
    :param verbose: verbosity level
    :return: a dict: for each station id (key), we get lists of dates, height, uncertainty data
    """
    # Build the GET params for the WFS request. We will use a CSV output, compact and easy to load into pandas
    sv_ids_str = "','".join(sv_names)
    params = {
        "service": "WFS",
        "request": "GetFeature",
        "version": "1.1.1",
        "typename": "hydroweb_rivers_ope_wse",
        "CQL_FILTER": f"name IN ('{sv_ids_str}') and start_time BETWEEN '{min_date.isoformat()}' and '{max_date.isoformat()}'",
        "outputFormat": "csv",
        "sortBy": "start_time"
    }
    encoded_params = urllib.parse.urlencode(params)

    req_url = f"{hydroweb_next_service_url}{encoded_params}"
    if verbose >= 1:
        print('URL {}'.format(req_url))
    # directly opening with read_csv doesn't seem to work, so using here the more robust traditional approach
    csv_data = requests.get(req_url).content
    df = pd.read_csv(io.StringIO(csv_data.decode('utf-8')), usecols=["name", "start_time", "wse", "wse_u"])
    # Parse dates as datetime. And since Pandas uses its own internal format, and the next function (and tests)
    # expects datetime, we convert it back to python datetime objects
    times = pd.to_datetime(df["start_time"]).dt.to_pydatetime()
    df["start_time"] = pd.Series(times, dtype="object")
    # Regroup the data by virtual station, then build the dicts as they were on previous function
    dg = df.groupby("name")
    values = {k: {"dates": t["start_time"].to_list(), "wsh": t["wse"].to_list(),
              "wsh_uncertainty": t["wse_u"].to_list(), } for k, t in dg}
    return values


def get_hysope_flowrates(sv_dict, min_date=None, max_date=None, verbose=1):

    # wsh_info = get_hysope_wsh(list(sv_dict.keys()), min_date=min_date, max_date=max_date, verbose=verbose)
    wsh_info = get_hysope_wsh_hydrowebnext(list(sv_dict.keys()), min_date=min_date, max_date=max_date, verbose=verbose)

    dico_out = dict()
    for sv_name, sv_info in sv_dict.items():
        dico_out[sv_name] = {'dates': wsh_info[sv_name]['dates']}
        flowrates_loc = np.ma.masked_invalid((np.array(wsh_info[sv_name]['wsh'], dtype=np.float64)-float(sv_info['coeff_Z0'])))
        flowrates_loc.mask[flowrates_loc<=0.] = True
        flowrates_loc = np.ma.masked_invalid(float(sv_info['coeff_A'])*(flowrates_loc**float(sv_info['coeff_B'])))
        if np.any(flowrates_loc.mask) and verbose >= 1:
            print('Invalid values in SV %s for dates :\n%s'%(sv_name, '\n'.join(['  ' + wsh_info[sv_name]['dates'][ii].strftime('%Y-%m-%dT%H:%M:%S') \
                for ii in range(len(wsh_info[sv_name]['dates'])) if flowrates_loc.mask[ii]])))
        dico_out[sv_name] = {'dates': [wsh_info[sv_name]['dates'][ii] for ii in range(len(wsh_info[sv_name]['dates'])) if not flowrates_loc.mask[ii]], \
            'flowrate': [flowrates_loc[ii] for ii in range(len(wsh_info[sv_name]['dates'])) if not flowrates_loc.mask[ii]]}
    return dico_out


def get_hysope_assimilation_data(sv_dict, min_date=None, max_date=None, verbose=1):
    flowrate_info = get_hysope_flowrates(sv_dict, min_date=min_date, max_date=max_date, verbose=verbose)
    meas_list = []
    for sv_name, sv_data in flowrate_info.items():
        for i0 in range(len(sv_data['dates'])):
            meas_list.append({'mesh_id': sv_dict[sv_name]['mesh_id'], 'sv_name': sv_name, 'date_data': sv_data['dates'][i0], \
                'value': sv_data['flowrate'][i0], 'uncertainty': sv_data['flowrate'][i0]*sv_dict[sv_name]['ratio_uncertainty']})
    return meas_list
