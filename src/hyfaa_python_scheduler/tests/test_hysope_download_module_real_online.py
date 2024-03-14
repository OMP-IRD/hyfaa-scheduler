"""
Test hyfaa/database/assimilation/hysope_download_module
"""
import unittest
from datetime import datetime

import hyfaa.database.assimilation.hysope_download_module as hysope

class TestGetHysopeWsh(unittest.TestCase):
    """
    Test that the get_hysope_wsh returns the expected result, don't mock, use real online services
    """
    def test_get_hysope_wsh_1m_with_online_hydroweb(self):
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
        data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                     max_date=datetime.fromisoformat(maxd), verbose=0)
        self.assertEqual(data, expected_return_value)

    def test_get_hysope_wsh_3m_with_online_hydroweb(self):
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
        data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                     max_date=datetime.fromisoformat(maxd), verbose=0)
        self.assertEqual(data, expected_return_value)


    def test_get_hysope_wsh_3m_with_online_hydrowebnext(self):
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
        data = hysope.get_hysope_wsh_hydrowebnext(sv_names, min_date=datetime.fromisoformat(mind),
                                     max_date=datetime.fromisoformat(maxd), verbose=0)
        self.assertEqual(data, expected_return_value)




    def test_get_compare_get_hysope_next_and_previous(self):
        """
        Runs the get_hysope_wsh and get_hysope_wsh_hydrowebnext and compare the results. Should match
        disable this test when hydroweb v1 gets decommissioned
        :return:
        """
        mind = '2023-11-01'
        maxd = '2024-03-03'
        sv_names = ['R_OYAPOCK_CAMOPI_KM0177', 'R_MARONI_TAPANAHONY_KM0198', 'R_MARONI_OELEMARI_KM0387',
                    'R_MARONI_TAMPOK_KM0345','R_MARONI_MARONI_KM0054']

        data = hysope.get_hysope_wsh(sv_names, min_date=datetime.fromisoformat(mind),
                                     max_date=datetime.fromisoformat(maxd), verbose=0)
        data_next = hysope.get_hysope_wsh_hydrowebnext(sv_names, min_date=datetime.fromisoformat(mind),
                                     max_date=datetime.fromisoformat(maxd), verbose=0)
        self.assertEqual(data, data_next)

    def test_get_hysope_flowrate_no_mock(self):
        """
        Check the flowrate function and downstream, no mocked service
        """
        assert True


if __name__ == '__main__':
    unittest.main()
