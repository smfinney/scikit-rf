# -*- coding: utf-8 -*-
import unittest
import os

from skrf.media import GCPW
from skrf.network import Network
from skrf.frequency import Frequency
import skrf as rf
from numpy.testing import assert_array_almost_equal, run_module_suite


class GCPWTestCase(unittest.TestCase):

    # Test parameters
    cmp_precision = 1 # best we can do for now

    freq = rf.Frequency(start=100, stop=10000, npoints=101, unit='MHz')


    # Test values calculated using qucstrans tool
    test_data = [

        #Default values
        {
            'line_args': {},
            'z0': 49.711,
            'er': 2.98,
        },

        # Sample geometry without metallization thickness
        {
            'line_args': {'w': 1.6e-3, 's': 3e-4, 'h': 1.524e-3, 'ep_r': 4.4},
            'z0': 50.106,
            'er': 2.84,
        },

         # Sample geometry with metallization thickness
        {
            'line_args': {'w': 1.6e-3, 's': 3e-4, 'h': 1.524e-3, 'ep_r': 4.4, 't': 1e-5},
            'z0': 49.617,
            'er': 2.80,
        },
    ]

    def setUp(self):

        # self.files_dir = os.path.join(
        #     os.path.dirname(os.path.abspath(__file__)),
        #     'qucs_prj'
        #     )        
        # fname = os.path.join(self.files_dir, 'gcpw.s2p')
        # self.qucs_ntwk = rf.Network(fname)

        self.networks = []

        for d in self.test_data:

            d['line'] = GCPW(frequency=self.freq, **d['line_args'])
            self.networks.append(d)

    def test_Z0(self):
        """
        Test modeled characteristic impedances against reference values
        """
        for n in self.networks:
            assert_array_almost_equal(n['line'].Z0, n['z0'], decimal=self.cmp_precision)
        
    def test_eps_eff(self):
        """
        Test modeled effective permittivity against reference values
        """
        for n in self.networks:
            assert_array_almost_equal(n['line'].ep_re, n['er'], decimal=self.cmp_precision)  
    

if __name__ == "__main__":
    # Launch all tests
    run_module_suite()