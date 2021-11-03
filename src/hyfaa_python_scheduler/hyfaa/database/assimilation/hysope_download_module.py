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

from hyfaa.utils.fileenv import fileenv

hydroweb_credentials = {
    'user': fileenv('HYDROWEB_USER'),
    'pwd' : fileenv('HYDROWEB_PASSWORD'),
}
hydroweb_baserequest = 'http://hydroweb.theia-land.fr/hydroweb/authdownload?products='


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


def get_hysope_flowrates(sv_dict, min_date=None, max_date=None, verbose=1):

    wsh_info = get_hysope_wsh(list(sv_dict.keys()), min_date=min_date, max_date=max_date, verbose=verbose)

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
