#!/usr/bin/env python
# -*- coding: utf-8 -*-
import geopandas
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
from eodag import EODataAccessGateway, SearchResult, setup_logging
from fnmatch import fnmatch
import logging
import numpy as np
import pandas as pd
from shapely.geometry import Point
import urllib
import io
import os
from shutil import rmtree

from hyfaa.utils.fileenv import fileenv

hydroweb_credentials = {
    'user': fileenv('HYDROWEB_USER'),
    'pwd': fileenv('HYDROWEB_PASSWORD'),
}
hydroweb_baserequest = 'http://hydroweb.theia-land.fr/hydroweb/authdownload?products='

hydroweb_next_service_url = "https://hydroweb.next.theia-land.fr/geoserver/HYDROWEB/ows?"

hydroweb_next_stac_url = "https://hydroweb.next.theia-land.fr/api/v1/rs-catalog/stac/"


def get_hysope_wsh(sv_names, min_date=None, max_date=None, verbose=1):
    print('using hydroweb username {}'.format(hydroweb_credentials['user']))
    dico_out = dict()
    tstart0 = datetime.utcnow()
    for sv_name in set(sv_names):
        tstart = datetime.utcnow()
        req_txt = hydroweb_baserequest + sv_name + '&format=json'
        if min_date is not None:
            req_txt += '&sdate=%s' % min_date.strftime('%Y-%m-%d')
        if max_date is not None:
            req_txt += '&edate=%s' % max_date.strftime('%Y-%m-%d')
        req_txt += '&user=%s&pwd=%s' % (hydroweb_credentials['user'], hydroweb_credentials['pwd'])
        if verbose >= 1:
            print('URL {}'.format(req_txt))
        data_loc = requests.get(req_txt).json()['data']
        datetimes_loc = [datetime.strptime(el['date'] + 'T' + el['time'], '%Y-%m-%dT%H:%M') for el in data_loc]
        values_loc = [float(el['orthometric_height_of_water_surface_at_reference_position']) if el[
                                                                                                    'orthometric_height_of_water_surface_at_reference_position'] is not None else None
                      for el in data_loc]
        uncertainties_loc = [float(el['associated_uncertainty']) if el['associated_uncertainty'] is not None else None
                             for el in data_loc]

        # select only non NAN values
        ids_select = [ii for ii in range(len(datetimes_loc)) if
                      (values_loc[ii] is not None) or (uncertainties_loc[ii] is not None)]
        datetimes_loc = [datetimes_loc[ii] for ii in ids_select]
        values_loc = [values_loc[ii] for ii in ids_select]
        uncertainties_loc = [uncertainties_loc[ii] for ii in ids_select]

        dico_out[sv_name] = {'dates': datetimes_loc, 'wsh': values_loc, 'wsh_uncertainty': uncertainties_loc}
        if verbose >= 1:
            print('Downloaded assimilation data from SV %s in %s' % (sv_name, datetime.utcnow() - tstart))
    if verbose >= 1:
        print('Downloaded all assimilation data in %s' % (datetime.utcnow() - tstart0))
    return dico_out


def get_hysope_wsh_hydrowebnext_gs(sv_names, min_date=None, max_date=None, verbose=1):
    """
    *** DEPRECATED ***: their geoserver services don't support large list of station identifiers, making this
    basicly useless
    Replaced by get_hysope_wsh_hydrowebnext


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


def get_name_from_hydrowebnext_filename(filename) -> str:
    """
    On hydrowebnext, the station name is no more easy to get. It is using an URN-like ID almost everywhere, except on
    the file's name, which is based on the old name, the ones we have in config/hyseop_svs.yaml config file,
    but slightly different
    :param filename:
    :return: station name (str)
    """
    # drop prefix and suffix
    filename_chunks = filename.split("_")
    filename_chunks = filename_chunks[1:-1]
    station_name = "_".join(filename_chunks)
    return station_name


def get_hysope_wsh_hydrowebnext(sv_dict, min_date=None, max_date=None, verbose=1):
    """
    Get water height from Hydroweb.next services. Should replace get_hysope_wsh function, that was using the old
    hydroweb API.

    Uses the official STAC API. But it is quite limited in functionalities, forcing us to work around those limitations,
    making this code much more complex than it should be...
    Uses the script sample provided by hydroweb.next platform (v 2.11+)

    It will expect that you provide the API Key as an environment variable
    named EODAG__HYDROWEB_NEXT__AUTH__CREDENTIALS__APIKEY

    :param sv_dict: list of virtual stations with attributes (interesting: coordinates)
    :param min_date: datetime object
    :param max_date: datetime object
    :param verbose: verbosity level
    :return: a dict, wsh_values: for each station id (key), we get lists of dates, height, uncertainty data
    """
    wsh_values = dict()  # will store the final results returned by this function

    # We compute the polygon enveloping the whole set of stations (points). We get is as WKT, which can be taken
    # as argument in the STAC search. We apply a buffer because if we don't some stations on the periphery
    # can get omitted
    pd_stations = geopandas.GeoSeries([Point(s["lon"], s["lat"]) for s in sv_dict.values()])
    stations_envelope_wkt = pd_stations.buffer(0.05).unary_union.envelope.wkt
    # stations_envelope_geojson = shapely.to_geojson(shapely.from_wkt(stations_envelope_wkt))
    # Using STAC client doesn't work. STAC implementation of hydrowebnext seems buggy, retrieving the search results
    # fails at some point making the whole process unusable. Last resort is using EODAG software

    # HYDROWEBNEXT_API_KEY = os.environ.get("HYDROWEBNEXT_API_KEY")
    # headers = {
    #     "X-API-Key": HYDROWEBNEXT_API_KEY,
    # }
    #
    # catalog = pystac_client.Client.open(hydroweb_next_stac_url, headers=headers)

    # results = catalog.search(
    #     collections="HYDROWEB_RIVERS_OPE",
    #     intersects=stations_envelope_geojson,
    #     #datetime=[min_date,  max_date]
    # )
    # for item in results.items():
    #     try:
    #         print(item.id)
    #     except Exception as e:
    #         pass

    setup_logging(0) # 0: nothing, 1: only progress bars, 2: INFO, 3: DEBUG
    # Add custom logger
    custom_logger = logging.getLogger("Hydroweb.next")
    custom_logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)-15s %(name)-32s [%(levelname)-8s] %(message)s')
    handler = logging.StreamHandler()  # Use a different handler if needed
    handler.setFormatter(formatter)
    custom_logger.addHandler(handler)

    custom_logger.info(f"Start retrieving data from hydroweb.next API")

    dag = EODataAccessGateway()

    # Set timeout to 30s
    # os.environ["EODAG__HYDROWEB_NEXT__SEARCH__TIMEOUT"] = "30"

    # Default search criteria when iterating over collection pages
    default_search_criteria = {
        "items_per_page": 2000,
    }

    # Set download directory
    tmpdir = os.environ.get("hyfaa_temp_dir", "/tmp")
    hysope_results_download_path = f"{tmpdir}/hysope_results_{int(datetime.now().timestamp())}"
    os.makedirs(hysope_results_download_path, exist_ok=True)

    query_args = {
        "productType": "HYDROWEB_RIVERS_OPE",
        "start": min_date.isoformat(),
        "end": max_date.isoformat(),
        "geom": stations_envelope_wkt,
    }
    query_args.update(default_search_criteria)

    # Run a paginated search
    search_results = ([])

    for i, page_results in enumerate(dag.search_iter_page(**query_args)):
        custom_logger.info(f"{len(page_results)} product(s) found on page {i + 1}")
        search_results.extend(page_results)

    custom_logger.info(f"Total products found : {len(search_results)}")

    # document the mapping between the old station ID and the new URN-like ID
    # And filter-out stations that we are not interested in (keep only those in sv_dict)
    ids_mappings = {}
    filtered_results = []
    for product in search_results:
        file_name = list(product.assets.data.keys())[0]
        station_name = get_name_from_hydrowebnext_filename(file_name)
        custom_logger.info(f"{station_name}      ({file_name})")
        ids_mappings[product.properties["id"]] = {
            "filename": file_name,
            "station_name": station_name
        }
        if station_name in list(sv_dict.keys()):
            filtered_results.append(product)

    custom_logger.info(f"Only {len(filtered_results)} products match the given list (hyfaa hysope config see config/hysope_svs.yaml file)")
    # Download all found products
    custom_logger.info(f"Downloading {len(filtered_results)} products...")
    downloaded_paths = dag.download_all(filtered_results, outputs_prefix=hysope_results_download_path)
    if downloaded_paths:
        distinct_values = list(set(downloaded_paths))
        # Check is distinct values length is equal to all results length
        if len(distinct_values) != len(filtered_results):
            custom_logger.warning(
                f"Distinct values length is not equal to all results length. {len(distinct_values)} != {len(filtered_results)}")
            custom_logger.warning(f"Some files have not been downloaded")
        else:
            custom_logger.info(f"All {len(filtered_results)} files have been successfully downloaded.")
    else:
        print(f"No files downloaded! Verify API-KEY and/or product search configuration.")

    # Process the downloaded datasets:
    # - filter by date
    # - store the results in wsh_values dict
    filename_pattern = "hydroprd_*_exp.txt"
    for path, subdirs, files in os.walk(hysope_results_download_path):
        for name in files:
            if fnmatch(name, filename_pattern):
                df = pd.read_csv(os.path.join(path, name), header=None, skiprows=45, sep=" ",
                                 names=["d", "t", "h", "h_u", "sep", "lon", "lat", "c7", "c8", "c9", "c10", "c11",
                                        "c12",
                                        "c13", "c14",
                                        "c15"])
                # Build a proper datetime column for filtering
                df["date"] = df[["d", "t"]].agg('T'.join, axis=1)
                df["date"] = pd.to_datetime(df['date'], format='%Y-%m-%dT%H:%M')
                # Filter data based on min-max date
                date_mask = (df['date'] >= min_date) & (df['date'] < max_date)
                filtered_df = df.loc[date_mask]

                wsh_values[get_name_from_hydrowebnext_filename(name)] = {
                    'dates': filtered_df["date"].tolist(),
                    'wsh': filtered_df["h"].tolist(),
                    'wsh_uncertainty': filtered_df["h_u"].tolist()
                }
    # Sort the keys in alphabetical order. Not necessary but easier to inspect/debug
    wsh_values_sorted = dict(sorted(wsh_values.items()))
    rmtree(hysope_results_download_path, ignore_errors=True)
    return wsh_values_sorted


def get_hysope_flowrates(sv_dict, min_date=None, max_date=None, verbose=1):
    # wsh_info = get_hysope_wsh(list(sv_dict.keys()), min_date=min_date, max_date=max_date, verbose=verbose)
    # wsh_info_gs = get_hysope_wsh_hydrowebnext_gs(list(sv_dict.keys()), min_date=min_date, max_date=max_date,
    #                       verbose=verbose)
    wsh_info = get_hysope_wsh_hydrowebnext(sv_dict, min_date=min_date, max_date=max_date, verbose=verbose)

    dico_out = dict()

    # Check for invalid values
    for sv_name, sv_info in sv_dict.items():
        dico_out[sv_name] = {'dates': wsh_info[sv_name]['dates']}
        flowrates_loc = np.ma.masked_invalid(
            (np.array(wsh_info[sv_name]['wsh'], dtype=np.float64) - float(sv_info['coeff_Z0'])))
        flowrates_loc.mask[flowrates_loc <= 0.] = True
        flowrates_loc = np.ma.masked_invalid(float(sv_info['coeff_A']) * (flowrates_loc ** float(sv_info['coeff_B'])))
        if np.any(flowrates_loc.mask) and verbose >= 1:
            print('Invalid values in SV %s for dates :\n%s' % (
                sv_name, '\n'.join(['  ' + wsh_info[sv_name]['dates'][ii].strftime('%Y-%m-%dT%H:%M:%S') \
                                    for ii in range(len(wsh_info[sv_name]['dates'])) if flowrates_loc.mask[ii]])))
        dico_out[sv_name] = {
            'dates': [wsh_info[sv_name]['dates'][ii] for ii in range(len(wsh_info[sv_name]['dates'])) if
                      not flowrates_loc.mask[ii]], \
            'flowrate': [flowrates_loc[ii] for ii in range(len(wsh_info[sv_name]['dates'])) if
                         not flowrates_loc.mask[ii]]}
    return dico_out


def get_hysope_assimilation_data(sv_dict, min_date=None, max_date=None, verbose=1):
    flowrate_info = get_hysope_flowrates(sv_dict, min_date=min_date, max_date=max_date, verbose=verbose)
    meas_list = []
    for sv_name, sv_data in flowrate_info.items():
        for i0 in range(len(sv_data['dates'])):
            meas_list.append(
                {'mesh_id': sv_dict[sv_name]['mesh_id'], 'sv_name': sv_name, 'date_data': sv_data['dates'][i0], \
                 'value': sv_data['flowrate'][i0],
                 'uncertainty': sv_data['flowrate'][i0] * sv_dict[sv_name]['ratio_uncertainty']})
    return meas_list
