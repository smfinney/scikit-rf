# -*- coding: utf-8 -*-
import unittest
import os

from skrf.media import GCPW
from skrf.network import Network
from skrf.frequency import Frequency
import skrf as rf
import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite


class GCPWTestCase(unittest.TestCase):

    # Test parameters

    freq_ext = rf.Frequency(start=100, stop=10000, npoints=101, unit='MHz')

    # Test values calculated using qucstrans tool
    freq_pts = rf.Frequency.from_f([0.01, 0.1, 1, 10], unit='GHz')
    test_data = [

        #Default values
        {
            'line_args': {},
            'z0': [49.7116, 49.7114, 49.6968, 48.8589],
            'er': [2.984, 2.984, 2.986, 3.089],
        },

        # Sample geometry without metallization thickness
        {
            'line_args': {'w': 1.6e-3, 's': 3e-4, 'h': 1.524e-3, 'ep_r': 4.4},
            'z0': [50.1058, 50.1056, 50.0962, 49.5328],
            'er': [2.836, 2.836, 2.837, 2.902],
        },

         # Sample geometry with metallization thickness
        {
            'line_args': {'w': 1.6e-3, 's': 3e-4, 'h': 1.524e-3, 'ep_r': 4.4, 't': 1e-5},
            'z0': [49.6170, 49.6169, 49.6072, 49.0290],
            'er': [2.795, 2.795, 2.796, 2.863],
        },
    ]

    def setUp(self):

        self.files_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'qucs_prj'
            )        
        fname = os.path.join(self.files_dir, 'gcpw.s2p')
        self.qucs_ntwk = rf.Network(fname)

        self.networks = []

        for d in self.test_data:

            d['line'] = GCPW(frequency=self.freq_pts, **d['line_args'])
            self.networks.append(d)

    def test_Z0(self):
        """
        Test modeled characteristic impedances against reference values
        """
        for n in self.networks:
            assert_array_almost_equal(n['line'].Z0, n['z0'], decimal=1)
        
    def test_eps_eff(self):
        """
        Test modeled effective permittivity against reference values
        """
        for n in self.networks:
            assert_array_almost_equal(n['line'].ep_re, n['er'], decimal=3)


    def test_s(self):
        """
        Test S-parameters against qucs model.
        """
        gcpw_model = GCPW(
            frequency=self.freq_ext,
            w=0.5e-3,
            s=0.1e-3,
            h=1e-3,
            ep_r=5,
            t=35e-6,
            rho=2.44e-8,
            )

        gcpw_ntwk = gcpw_model.thru(z0=50) ** gcpw_model.line(1, unit='m') ** gcpw_model.thru(z0=50)

        #print(gcpw_ntwk.s_mag[:,0,0])
        #print(self.qucs_ntwk.s_mag[:,0,0])

        # S11
        assert_array_almost_equal(gcpw_ntwk.s_mag[:,0,0], self.qucs_ntwk.s_mag[:,0,0], decimal=1)

        # S21
        #assert_array_almost_equal(gcpw_ntwk.s_mag[:,1,0], self.qucs_ntwk.s_mag[:,1,0], decimal=1)
    

if __name__ == "__main__":
    # Launch all tests
    run_module_suite()