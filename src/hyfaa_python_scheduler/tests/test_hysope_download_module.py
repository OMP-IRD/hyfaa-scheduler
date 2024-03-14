"""
Test hyfaa/database/assimilation/hysope_download_module
"""
import unittest
from datetime import datetime
from unittest.mock import patch

import requests_mock

import hyfaa.database.assimilation.hysope_download_module as hysope


class TestGetHysopeWsh(unittest.TestCase):
    """
    Test that the get_hysope_wsh returns the expected result
    """
    def test_get_hysope_wsh_1m_with_mock_hydroweb(self):
        with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
            mind = '2024-01-01'
            maxd = '2024-01-31'
            sv_names = ['R_APPROUAGUE_APPROUAGUE_KM0041', 'R_APPROUAGUE_APPROUAGUE_KM0047', 'R_MARONI_MARONI_KM0054']
            expected_return_value = {
                'R_MARONI_MARONI_KM0054': {'dates': [datetime(2024, 1, 29, 13, 40)], 'wsh': [1.38],
                                           'wsh_uncertainty': [0.02]},
                'R_APPROUAGUE_APPROUAGUE_KM0047': {'dates': [datetime(2024, 1, 19, 1, 26)], 'wsh': [0.84],
                                                   'wsh_uncertainty': [0.03]},
                'R_APPROUAGUE_APPROUAGUE_KM0041': {'dates': [datetime(2024, 1, 24, 13, 31)], 'wsh': [0.72],
                                                   'wsh_uncertainty': [0.01]}}
            with requests_mock.Mocker() as rm:
                for s in sv_names:
                    req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={mind}&edate={maxd}&user=&pwd='
                    with open(
                            f'data/hydroprd_{s}_from_{mind.replace("-", "")}_to_{maxd.replace("-", "")}.json') as f_in:
                        mock_data = f_in.read()
                        rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)

    def test_get_hysope_wsh_3m_with_mock_hydroweb(self):
        with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
            mind = '2023-11-01'
            maxd = '2024-01-31'
            sv_names = ['R_OYAPOCK_CAMOPI_KM0177', 'R_MARONI_TAPANAHONY_KM0198', 'R_MARONI_OELEMARI_KM0387',
                        'R_MARONI_TAMPOK_KM0345']
            expected_return_value = {
                'R_OYAPOCK_CAMOPI_KM0177': {
                    'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                              datetime(2023, 12, 28, 13, 32), datetime(2024, 1, 24, 13, 32)],
                    'wsh': [54.77, 55.02, 55.45, 56.16],
                    'wsh_uncertainty': [0.07, 0.01, 0.04, 0.01]
                },
                'R_MARONI_TAMPOK_KM0345': {
                    'dates': [datetime(2023, 11, 18, 1, 34),
                              datetime(2023, 12, 15, 1, 34),
                              datetime(2024, 1, 11, 1, 33)],
                    'wsh': [94.56, 95.44, 95.7],
                    'wsh_uncertainty': [0.03, 0.01, 0.02]
                },
                'R_MARONI_OELEMARI_KM0387': {
                    'dates': [datetime(2023, 11, 23, 13, 39),
                              datetime(2023, 12, 20, 13, 39),
                              datetime(2024, 1, 16, 13, 39)],
                    'wsh': [122.36, 122.89, 123.08],
                    'wsh_uncertainty': [0.22, 0.22, 0.04]
                },
                'R_MARONI_TAPANAHONY_KM0198': {
                    'dates': [datetime(2023, 11, 4, 1, 35),
                              datetime(2023, 12, 1, 1, 35),
                              datetime(2023, 12, 28, 1, 35),
                              datetime(2024, 1, 24, 1, 35)],
                    'wsh': [48.96, 49.0, 49.15, 48.59],
                    'wsh_uncertainty': [0.06, 0.02, 0.08, 0.11]}
            }
            with requests_mock.Mocker() as rm:
                for s in sv_names:
                    req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={mind}&edate={maxd}&user=&pwd='
                    with open(
                            f'data/hydroprd_{s}_from_{mind.replace("-", "")}_to_{maxd.replace("-", "")}.json') as f_in:
                        mock_data = f_in.read()
                        rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)

    def test_get_hysope_wsh_3m_with_mock_hydroweb(self):
        with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
            mind = '2023-11-01'
            maxd = '2024-01-31'
            sv_names = ['R_OYAPOCK_CAMOPI_KM0177', 'R_MARONI_TAPANAHONY_KM0198', 'R_MARONI_OELEMARI_KM0387',
                        'R_MARONI_TAMPOK_KM0345']
            expected_return_value = {
                'R_OYAPOCK_CAMOPI_KM0177': {
                    'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                              datetime(2023, 12, 28, 13, 32), datetime(2024, 1, 24, 13, 32)],
                    'wsh': [54.77, 55.02, 55.45, 56.16],
                    'wsh_uncertainty': [0.07, 0.01, 0.04, 0.01]
                },
                'R_MARONI_TAMPOK_KM0345': {
                    'dates': [datetime(2023, 11, 18, 1, 34),
                              datetime(2023, 12, 15, 1, 34),
                              datetime(2024, 1, 11, 1, 33)],
                    'wsh': [94.56, 95.44, 95.7],
                    'wsh_uncertainty': [0.03, 0.01, 0.02]
                },
                'R_MARONI_OELEMARI_KM0387': {
                    'dates': [datetime(2023, 11, 23, 13, 39),
                              datetime(2023, 12, 20, 13, 39),
                              datetime(2024, 1, 16, 13, 39)],
                    'wsh': [122.36, 122.89, 123.08],
                    'wsh_uncertainty': [0.22, 0.22, 0.04]
                },
                'R_MARONI_TAPANAHONY_KM0198': {
                    'dates': [datetime(2023, 11, 4, 1, 35),
                              datetime(2023, 12, 1, 1, 35),
                              datetime(2023, 12, 28, 1, 35),
                              datetime(2024, 1, 24, 1, 35)],
                    'wsh': [48.96, 49.0, 49.15, 48.59],
                    'wsh_uncertainty': [0.06, 0.02, 0.08, 0.11]}
            }
            with requests_mock.Mocker() as rm:
                for s in sv_names:
                    req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={mind}&edate={maxd}&user=&pwd='
                    with open(
                            f'data/hydroprd_{s}_from_{mind.replace("-", "")}_to_{maxd.replace("-", "")}.json') as f_in:
                        mock_data = f_in.read()
                        rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)


class TestGetHysopeWshNext(unittest.TestCase):
    """
    Test that the get_hysope_wsh_hydrowebnext returns the expected result
    Run the same tests as for the previous version (see above)
    """
    def test_get_hysope_wsh_next_1m_with_mock_hydroweb_next(self):
        mind = '2024-01-01'
        maxd = '2024-01-31'
        sv_names = ['R_APPROUAGUE_APPROUAGUE_KM0041', 'R_APPROUAGUE_APPROUAGUE_KM0047', 'R_MARONI_MARONI_KM0054']
        expected_return_value = {
            'R_MARONI_MARONI_KM0054': {'dates': [datetime(2024, 1, 29, 13, 40)], 'wsh': [1.38],
                                       'wsh_uncertainty': [0.02]},
            'R_APPROUAGUE_APPROUAGUE_KM0047': {'dates': [datetime(2024, 1, 19, 1, 26)], 'wsh': [0.84],
                                               'wsh_uncertainty': [0.03]},
            'R_APPROUAGUE_APPROUAGUE_KM0041': {'dates': [datetime(2024, 1, 24, 13, 31)], 'wsh': [0.72],
                                               'wsh_uncertainty': [0.01]}}
        with requests_mock.Mocker() as rm:
            req_url = f'https://hydroweb.next.theia-land.fr/geoserver/HYDROWEB/ows?service=WFS&request=GetFeature&version=1.1.1&typename=hydroweb_rivers_wse&CQL_FILTER=name+IN+%28%27R_APPROUAGUE_APPROUAGUE_KM0041%27%2C%27R_APPROUAGUE_APPROUAGUE_KM0047%27%2C%27R_MARONI_MARONI_KM0054%27%29+and+start_time+BETWEEN+%272024-01-01T00%3A00%3A00%27+and+%272024-01-31T00%3A00%3A00%27&outputFormat=csv&sortBy=start_time'
            with open(f'data/hydroweb_rivers_wse_2024-01-01-to-2024-01-31.csv') as f_in:
                mock_data = f_in.read()
                rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh_hydrowebnext(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)

    def test_get_hysope_wsh_3m_with_mock_hydroweb(self):
        with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
            mind = '2023-11-01'
            maxd = '2024-01-31'
            sv_names = ['R_OYAPOCK_CAMOPI_KM0177', 'R_MARONI_TAPANAHONY_KM0198', 'R_MARONI_OELEMARI_KM0387',
                        'R_MARONI_TAMPOK_KM0345']
            expected_return_value = {
                'R_OYAPOCK_CAMOPI_KM0177': {
                    'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                              datetime(2023, 12, 28, 13, 32), datetime(2024, 1, 24, 13, 32)],
                    'wsh': [54.77, 55.02, 55.45, 56.16],
                    'wsh_uncertainty': [0.07, 0.01, 0.04, 0.01]
                },
                'R_MARONI_TAMPOK_KM0345': {
                    'dates': [datetime(2023, 11, 18, 1, 34),
                              datetime(2023, 12, 15, 1, 34),
                              datetime(2024, 1, 11, 1, 33)],
                    'wsh': [94.56, 95.44, 95.7],
                    'wsh_uncertainty': [0.03, 0.01, 0.02]
                },
                'R_MARONI_OELEMARI_KM0387': {
                    'dates': [datetime(2023, 11, 23, 13, 39),
                              datetime(2023, 12, 20, 13, 39),
                              datetime(2024, 1, 16, 13, 39)],
                    'wsh': [122.36, 122.89, 123.08],
                    'wsh_uncertainty': [0.22, 0.22, 0.04]
                },
                'R_MARONI_TAPANAHONY_KM0198': {
                    'dates': [datetime(2023, 11, 4, 1, 35),
                              datetime(2023, 12, 1, 1, 35),
                              datetime(2023, 12, 28, 1, 35),
                              datetime(2024, 1, 24, 1, 35)],
                    'wsh': [48.96, 49.0, 49.15, 48.59],
                    'wsh_uncertainty': [0.06, 0.02, 0.08, 0.11]}
            }
            with requests_mock.Mocker() as rm:
                for s in sv_names:
                    req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={mind}&edate={maxd}&user=&pwd='
                    with open(
                            f'data/hydroprd_{s}_from_{mind.replace("-", "")}_to_{maxd.replace("-", "")}.json') as f_in:
                        mock_data = f_in.read()
                        rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)

    def test_get_hysope_wsh_3m_with_mock_hydroweb(self):
        with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
            mind = '2023-11-01'
            maxd = '2024-01-31'
            sv_names = ['R_OYAPOCK_CAMOPI_KM0177', 'R_MARONI_TAPANAHONY_KM0198', 'R_MARONI_OELEMARI_KM0387',
                        'R_MARONI_TAMPOK_KM0345']
            expected_return_value = {
                'R_OYAPOCK_CAMOPI_KM0177': {
                    'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                              datetime(2023, 12, 28, 13, 32), datetime(2024, 1, 24, 13, 32)],
                    'wsh': [54.77, 55.02, 55.45, 56.16],
                    'wsh_uncertainty': [0.07, 0.01, 0.04, 0.01]
                },
                'R_MARONI_TAMPOK_KM0345': {
                    'dates': [datetime(2023, 11, 18, 1, 34),
                              datetime(2023, 12, 15, 1, 34),
                              datetime(2024, 1, 11, 1, 33)],
                    'wsh': [94.56, 95.44, 95.7],
                    'wsh_uncertainty': [0.03, 0.01, 0.02]
                },
                'R_MARONI_OELEMARI_KM0387': {
                    'dates': [datetime(2023, 11, 23, 13, 39),
                              datetime(2023, 12, 20, 13, 39),
                              datetime(2024, 1, 16, 13, 39)],
                    'wsh': [122.36, 122.89, 123.08],
                    'wsh_uncertainty': [0.22, 0.22, 0.04]
                },
                'R_MARONI_TAPANAHONY_KM0198': {
                    'dates': [datetime(2023, 11, 4, 1, 35),
                              datetime(2023, 12, 1, 1, 35),
                              datetime(2023, 12, 28, 1, 35),
                              datetime(2024, 1, 24, 1, 35)],
                    'wsh': [48.96, 49.0, 49.15, 48.59],
                    'wsh_uncertainty': [0.06, 0.02, 0.08, 0.11]}
            }
            with requests_mock.Mocker() as rm:
                for s in sv_names:
                    req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={mind}&edate={maxd}&user=&pwd='
                    with open(
                            f'data/hydroprd_{s}_from_{mind.replace("-", "")}_to_{maxd.replace("-", "")}.json') as f_in:
                        mock_data = f_in.read()
                        rm.get(req_url, text=mock_data)

                data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                             max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)



class TestGetHysopeFlowrates(unittest.TestCase):
    """
    Test get_hysope_flowrates
    """

    def test_get_hysope_flowrate_mock_wsh(self):
        """
        Mock the get_hysope_wsh function, check the flowrate function itself
        """
        mind = '2023-11-01'
        maxd = '2024-01-31'
        sv_dict = {
            'R_MARONI_OELEMARI_KM0387': {'coeff_A': 92.88, 'coeff_B': 1.077, 'coeff_Z0': 122.496, 'lat': 3.194,
                                         'lon': -54.286, 'mesh_id': 4153, 'ratio_uncertainty': 0.1398},
            'R_MARONI_TAMPOK_KM0345': {'coeff_A': 4.927, 'coeff_B': 2.325, 'coeff_Z0': 92.348, 'lat': 3.398,
                                       'lon': -53.894, 'mesh_id': 4264, 'ratio_uncertainty': 0.0753},
            'R_MARONI_TAPANAHONY_KM0198': {'coeff_A': 398.039, 'coeff_B': 1.851, 'coeff_Z0': 48.17, 'lat': 4.229,
                                           'lon': -54.546, 'mesh_id': 4411, 'ratio_uncertainty': 0.1886},
            'R_OYAPOCK_CAMOPI_KM0177': {'coeff_A': 43.165, 'coeff_B': 1.578, 'coeff_Z0': 53.93, 'lat': 3.208,
                                        'lon': -52.412, 'mesh_id': 4224, 'ratio_uncertainty': 0.0388},
        }
        mock_wsh_value = {
            'R_OYAPOCK_CAMOPI_KM0177': {
                'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                          datetime(2023, 12, 28, 13, 32), datetime(2024, 1, 24, 13, 32)],
                'wsh': [54.77, 55.02, 55.45, 56.16],
                'wsh_uncertainty': [0.07, 0.01, 0.04, 0.01]
            },
            'R_MARONI_TAMPOK_KM0345': {
                'dates': [datetime(2023, 11, 18, 1, 34),
                          datetime(2023, 12, 15, 1, 34),
                          datetime(2024, 1, 11, 1, 33)],
                'wsh': [94.56, 95.44, 95.7],
                'wsh_uncertainty': [0.03, 0.01, 0.02]
            },
            'R_MARONI_OELEMARI_KM0387': {
                'dates': [datetime(2023, 11, 23, 13, 39),
                          datetime(2023, 12, 20, 13, 39),
                          datetime(2024, 1, 16, 13, 39)],
                'wsh': [122.36, 122.89, 123.08],
                'wsh_uncertainty': [0.22, 0.22, 0.04]
            },
            'R_MARONI_TAPANAHONY_KM0198': {
                'dates': [datetime(2023, 11, 4, 1, 35),
                          datetime(2023, 12, 1, 1, 35),
                          datetime(2023, 12, 28, 1, 35),
                          datetime(2024, 1, 24, 1, 35)],
                'wsh': [48.96, 49.0, 49.15, 48.59],
                'wsh_uncertainty': [0.06, 0.02, 0.08, 0.11]}
        }
        expected_return_value = {
            'R_MARONI_OELEMARI_KM0387': {
                'dates': [datetime(2023, 12, 20, 13, 39), datetime(2024, 1, 16, 13, 39)],
                'flowrate': [34.06211709933239, 52.04138463564574]
            },
            'R_MARONI_TAMPOK_KM0345': {
                'dates': [datetime(2023, 11, 18, 1, 34), datetime(2023, 12, 15, 1, 34),
                          datetime(2024, 1, 11, 1, 33)],
                'flowrate': [31.203834973260378, 67.98128437390076, 82.01898404005698]
            },
            'R_MARONI_TAPANAHONY_KM0198': {
                'dates': [datetime(2023, 11, 4, 1, 35), datetime(2023, 12, 1, 1, 35),
                          datetime(2023, 12, 28, 1, 35), datetime(2024, 1, 24, 1, 35)],
                'flowrate': [257.2962001697569, 281.9286264479068, 383.42911979029026, 79.90244392587749]
            },
            'R_OYAPOCK_CAMOPI_KM0177': {
                'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                          datetime(2023, 12, 28, 13, 32),
                          datetime(2024, 1, 24, 13, 32)],
                'flowrate': [32.782680315706365, 49.452785512250415, 83.57587401416401,
                             153.02303774625886]
            }
        }
        with patch('hyfaa.database.assimilation.hysope_download_module.get_hysope_wsh', return_value=mock_wsh_value):
            with patch('hyfaa.database.assimilation.hysope_download_module.get_hysope_wsh_hydrowebnext', return_value=mock_wsh_value):
                data = hysope.get_hysope_flowrates(sv_dict, min_date=datetime.fromisoformat(mind),
                                                   max_date=datetime.fromisoformat(maxd), verbose=0)
                self.assertEqual(data, expected_return_value)


class TestGetHysopeAssimilationData(unittest.TestCase):
    """
    Test get_hysope_assimilation_data
    """

    @classmethod
    def setUpClass(cls):
        cls.mind = '2023-11-01'
        cls.maxd = '2024-01-31'
        cls.sv_dict = {
            'R_MARONI_OELEMARI_KM0387': {'coeff_A': 92.88, 'coeff_B': 1.077, 'coeff_Z0': 122.496, 'lat': 3.194,
                                         'lon': -54.286, 'mesh_id': 4153, 'ratio_uncertainty': 0.1398},
            'R_MARONI_TAMPOK_KM0345': {'coeff_A': 4.927, 'coeff_B': 2.325, 'coeff_Z0': 92.348, 'lat': 3.398,
                                       'lon': -53.894, 'mesh_id': 4264, 'ratio_uncertainty': 0.0753},
            'R_MARONI_TAPANAHONY_KM0198': {'coeff_A': 398.039, 'coeff_B': 1.851, 'coeff_Z0': 48.17, 'lat': 4.229,
                                           'lon': -54.546, 'mesh_id': 4411, 'ratio_uncertainty': 0.1886},
            'R_OYAPOCK_CAMOPI_KM0177': {'coeff_A': 43.165, 'coeff_B': 1.578, 'coeff_Z0': 53.93, 'lat': 3.208,
                                        'lon': -52.412, 'mesh_id': 4224, 'ratio_uncertainty': 0.0388},
        }

        cls.expected_return_value = [
            {'mesh_id': 4153, 'sv_name': 'R_MARONI_OELEMARI_KM0387',
             'date_data': datetime(2023, 12, 20, 13, 39), 'value': 34.06211709933239,
             'uncertainty': 4.761883970486668},
            {'mesh_id': 4153, 'sv_name': 'R_MARONI_OELEMARI_KM0387',
             'date_data': datetime(2024, 1, 16, 13, 39), 'value': 52.04138463564574,
             'uncertainty': 7.275385572063275},
            {'mesh_id': 4264, 'sv_name': 'R_MARONI_TAMPOK_KM0345',
             'date_data': datetime(2023, 11, 18, 1, 34), 'value': 31.203834973260378,
             'uncertainty': 2.3496487734865066},
            {'mesh_id': 4264, 'sv_name': 'R_MARONI_TAMPOK_KM0345',
             'date_data': datetime(2023, 12, 15, 1, 34), 'value': 67.98128437390076,
             'uncertainty': 5.118990713354727},
            {'mesh_id': 4264, 'sv_name': 'R_MARONI_TAMPOK_KM0345',
             'date_data': datetime(2024, 1, 11, 1, 33), 'value': 82.01898404005698,
             'uncertainty': 6.176029498216291},
            {'mesh_id': 4411, 'sv_name': 'R_MARONI_TAPANAHONY_KM0198',
             'date_data': datetime(2023, 11, 4, 1, 35), 'value': 257.2962001697569,
             'uncertainty': 48.52606335201615},
            {'mesh_id': 4411, 'sv_name': 'R_MARONI_TAPANAHONY_KM0198',
             'date_data': datetime(2023, 12, 1, 1, 35), 'value': 281.9286264479068,
             'uncertainty': 53.17173894807522},
            {'mesh_id': 4411, 'sv_name': 'R_MARONI_TAPANAHONY_KM0198',
             'date_data': datetime(2023, 12, 28, 1, 35), 'value': 383.42911979029026,
             'uncertainty': 72.31473199244874},
            {'mesh_id': 4411, 'sv_name': 'R_MARONI_TAPANAHONY_KM0198',
             'date_data': datetime(2024, 1, 24, 1, 35), 'value': 79.90244392587749,
             'uncertainty': 15.069600924420493},
            {'mesh_id': 4224, 'sv_name': 'R_OYAPOCK_CAMOPI_KM0177',
             'date_data': datetime(2023, 11, 4, 13, 32), 'value': 32.782680315706365,
             'uncertainty': 1.271967996249407},
            {'mesh_id': 4224, 'sv_name': 'R_OYAPOCK_CAMOPI_KM0177',
             'date_data': datetime(2023, 12, 1, 13, 32), 'value': 49.452785512250415,
             'uncertainty': 1.9187680778753162},
            {'mesh_id': 4224, 'sv_name': 'R_OYAPOCK_CAMOPI_KM0177',
             'date_data': datetime(2023, 12, 28, 13, 32), 'value': 83.57587401416401,
             'uncertainty': 3.242743911749564},
            {'mesh_id': 4224, 'sv_name': 'R_OYAPOCK_CAMOPI_KM0177',
             'date_data': datetime(2024, 1, 24, 13, 32), 'value': 153.02303774625886,
             'uncertainty': 5.937293864554844}
        ]

    def test_get_hysope_assimilation_data_mock_flowrate(self):
        """
        Mock the flowrate function, check the get_hysope_assimilation_data function itself
        """
        mock_flowrate_value = {
            'R_MARONI_OELEMARI_KM0387': {
                'dates': [datetime(2023, 12, 20, 13, 39), datetime(2024, 1, 16, 13, 39)],
                'flowrate': [34.06211709933239, 52.04138463564574]
            },
            'R_MARONI_TAMPOK_KM0345': {
                'dates': [datetime(2023, 11, 18, 1, 34), datetime(2023, 12, 15, 1, 34),
                          datetime(2024, 1, 11, 1, 33)],
                'flowrate': [31.203834973260378, 67.98128437390076, 82.01898404005698]
            },
            'R_MARONI_TAPANAHONY_KM0198': {
                'dates': [datetime(2023, 11, 4, 1, 35), datetime(2023, 12, 1, 1, 35),
                          datetime(2023, 12, 28, 1, 35), datetime(2024, 1, 24, 1, 35)],
                'flowrate': [257.2962001697569, 281.9286264479068, 383.42911979029026, 79.90244392587749]
            },
            'R_OYAPOCK_CAMOPI_KM0177': {
                'dates': [datetime(2023, 11, 4, 13, 32), datetime(2023, 12, 1, 13, 32),
                          datetime(2023, 12, 28, 13, 32),
                          datetime(2024, 1, 24, 13, 32)],
                'flowrate': [32.782680315706365, 49.452785512250415, 83.57587401416401,
                             153.02303774625886]
            }
        }
        with patch('hyfaa.database.assimilation.hysope_download_module.get_hysope_flowrates',
                   return_value=mock_flowrate_value):
            data = hysope.get_hysope_assimilation_data(self.sv_dict, min_date=datetime.fromisoformat(self.mind),
                                                       max_date=datetime.fromisoformat(self.maxd), verbose=0)
            self.assertEqual(data, self.expected_return_value)


    # This one will need refactoring if it were to rely on the new hydroweb.next service
    # def test_get_hysope_assimilation_data_mock_hydroweb_service(self):
    #     """
    #     Mock the hydroweb online service (theia), check the full processing chain from the module
    #     """
    #     sv_names = self.sv_dict.keys()
    #     with patch.dict(hysope.hydroweb_credentials, {'user': '', 'pwd': ''}, clear=True):
    #         with requests_mock.Mocker() as rm:
    #             for s in sv_names:
    #                 req_url = f'http://hydroweb.theia-land.fr/hydroweb/authdownload?products={s}&format=json&sdate={self.mind}&edate={self.maxd}&user=&pwd='
    #                 with open(
    #                         f'data/hydroprd_{s}_from_{self.mind.replace("-", "")}_to_{self.maxd.replace("-", "")}.json') as f_in:
    #                     mock_data = f_in.read()
    #                     rm.get(req_url, text=mock_data)
    #             data = hysope.get_hysope_assimilation_data(self.sv_dict, min_date=datetime.fromisoformat(self.mind),
    #                                                        max_date=datetime.fromisoformat(self.maxd), verbose=0)
    #             self.assertEqual(data, self.expected_return_value)


if __name__ == '__main__':
    unittest.main()
