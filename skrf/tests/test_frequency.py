import os
import unittest
import warnings

import numpy as npy

import skrf as rf
from skrf.frequency import InvalidFrequencyWarning


class FrequencyTestCase(unittest.TestCase):
    '''

    '''
    def setUp(self):
        '''
        '''
        self.test_dir = os.path.dirname(os.path.abspath(__file__))+'/'

    def test_create_linear_sweep(self):
        freq = rf.Frequency(1,10,10,'ghz')
        self.assertTrue((freq.f == npy.linspace(1,10,10)*1e9).all())
        self.assertTrue((freq.f_scaled ==npy.linspace(1,10,10)).all())
        self.assertTrue((freq.sweep_type == 'lin'))

    def test_create_log_sweep(self):
        freq = rf.Frequency(1,10,10,'ghz', sweep_type='log')
        #Check end points
        self.assertTrue((freq.f[0] == 1e9))
        self.assertTrue((freq.f[-1] == 10e9))
        spacing = [freq.f[i+1]/freq.f[i] for i in range(len(freq.f)-1)]
        #Check that frequency is increasing
        self.assertTrue(all(s > 1 for s in spacing))
        #Check that ratio of adjacent frequency points is identical
        self.assertTrue(all(abs(spacing[i] - spacing[0]) < 1e-10 for i in range(len(spacing))))
        self.assertTrue((freq.sweep_type == 'log'))

    def test_create_rando_sweep(self):
        f = npy.array([1,5,200])
        freq = rf.Frequency.from_f(f,unit='khz')
        self.assertTrue((freq.f ==f*1e3).all())
        self.assertTrue((freq.f_scaled== f).all())
        self.assertTrue((freq.sweep_type == 'unknown'))
        
        with self.assertRaises(ValueError):
            freq.npoints = 10

    def test_rando_sweep_from_touchstone(self):
        '''
        this also tests the ability to read a touchstone file.
        '''
        rando_sweep_ntwk = rf.Network(os.path.join(self.test_dir, 'ntwk_arbitrary_frequency.s2p'))
        self.assertTrue((rando_sweep_ntwk.f == \
            npy.array([1,4,10,20])).all())
        self.assertTrue((rando_sweep_ntwk.frequency.sweep_type == 'unknown'))

    def test_slicer(self):
        a = rf.Frequency.from_f([1,2,4,5,6])

        b = a['2-5ghz']
        tinyfloat = 1e-12
        self.assertTrue((abs(b.f - [2e9,4e9,5e9]) < tinyfloat).all())

    def test_frequency_check(self):
        with self.assertWarns(InvalidFrequencyWarning):
            freq = rf.Frequency.from_f([2,1])

        with self.assertWarns(InvalidFrequencyWarning):
            freq = rf.Frequency.from_f([1,2,2])

        with warnings.catch_warnings(record=True) as warns:
            freq = rf.Frequency.from_f([1,2,2], unit="Hz")

            w = [w for w in warns if issubclass(w.category, InvalidFrequencyWarning)]
            if w:
                inv = freq.drop_non_monotonic_increasing()
                self.assertListEqual(inv, [2])
            
            else:
                self.fail("Warning is not triggered")
            
            self.assertTrue(npy.allclose(freq.f, [1,2]))




suite = unittest.TestLoader().loadTestsFromTestCase(FrequencyTestCase)
unittest.TextTestRunner(verbosity=2).run(suite)
