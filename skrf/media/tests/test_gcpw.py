# -*- coding: utf-8 -*-
import unittest
import os
import numpy as npy

from skrf.media import GCPW
from skrf.network import Network
from skrf.frequency import Frequency
import skrf as rf
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_allclose, run_module_suite


class GCPWTestCase(unittest.TestCase):

    def setUp(self):
        self.files_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'qucs_prj'
            )        
        fname = os.path.join(self.files_dir, 'gcpw.s2p')
        self.qucs_ntwk = rf.Network(fname)

        # create frequency
        self.freq = rf.Frequency(start=1, stop=20, npoints=21, unit='GHz')

        # create examples
        # # infinite dielectric substrate, infinitely thin metal
        # self.cpw1 = CPW(frequency=self.freq, w=40e-6, s=20e-6, ep_r=3)
        # # infinite GaAs substrate, infinitely thin metal
        # self.cpw2 = CPW(frequency=self.freq, w=75e-6, s=50e-6, ep_r=12.9)
        # # infinite GaAs substrate, finite metal thickess
        # # TODO: not used yet
        # self.cpw3 = CPW(frequency=self.freq, w=75e-6, s=50e-6, ep_r=12.9, t=1e-6)

    def test_qucs_network(self):
        """
        Test against the Qucs project results.

        Based on GaAs example from CPW test cases.

        Line properties:

        Width: 75 um (75e-6 m)
        Spacing: 50 um (50e-6 m)
        Substrate thickness: 1 mm (1e-3 m)
        Er: 12.9 (GaAs)

        Properties in non-ideal case:

        Metalization thickness: 0.01 um (1e-8 m)
        Rho: 2.44e-8 ohm-m (gold)
        Loss tangent: 2e-4
        

        TODO : finalize
        """
        # Create an equivalent test case without losses and metalization thickness
        gcpw = GCPW(
            frequency=self.freq, 
            w=75e-6, 
            s=50e-6, 
            h=1e-3, 
            ep_r=12.9,
            #tan_delta=2e-4,
            #t=1e-8, 
            #rho=2.44e-8,
            )

        ntw = gcpw.thru(z0=50)**gcpw.line(d=1, unit='m')**gcpw.thru(z0=50)
        #print(ntw.s)
        #print(50 * (1 + ntw.s11.s) / (1 - ntw.s11.s))

        #assert_array_equal(ntw.frequency.f, self.qucs_ntwk.frequency.f)

        assert_array_almost_equal(ntw.s, self.qucs_ntwk.s, decimal=3, verbose=True)
        #self.qucs_ntwk.plot_s_db()
        #ntw.plot_s_db()

    # def test_Z0(self):
    #     """
    #     Test the CPW Characteristic Impedances
    #     """
    #     assert_array_almost_equal(self.cpw1.Z0, 85.25, decimal=3)
        
    # def test_eps_eff(self):
    #     """
    #     Test the effective permittivity of CPW
    #     """
    #     assert_array_almost_equal(self.cpw1.ep_re, 2.00, decimal=3)
    #     assert_array_almost_equal(self.cpw2.ep_re, 6.95, decimal=3)        
        
    # def test_Z0_vs_f(self):
    #     """
    #     Test the CPW Characteristic Impedance vs frequency. 
        
    #     Reference data comes from Qucs Documentation (Fig 12.2)
    #     """        
    #     w_over_s_qucs, Z0_qucs = npy.loadtxt(
    #         os.path.join(self.files_dir, 'cpw_qucs_ep_r9dot5.csv'), 
    #         delimiter=';', unpack=True)
               
    #     w = 1
    #     Z0 = []
    #     for w_o_s in w_over_s_qucs:
    #         _cpw = CPW(frequency=self.freq[0], w=w, s=w/w_o_s, ep_r=9.5)
    #         Z0.append(_cpw.Z0[0].real)
            
    #     # all to a 3% relative difference
    #     # this is quite a large discrepency, but I extracted the ref values from the plot
    #     # one could do better eventually by extracting values from Qucs directly
    #     rel_diff = (Z0_qucs-npy.array(Z0))/Z0_qucs
    #     assert_allclose(rel_diff  - 3/100, 0, atol=0.1)

if __name__ == "__main__":
    # Launch all tests
    run_module_suite()