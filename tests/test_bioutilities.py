from unittest import TestCase
import numpy as np
from haystack.bioutilities import smooth

class Test(TestCase):
    def test_smooth(self):
        np.random.seed(33)
        raw_values = np.random.binomial(n=20, p=0.5, size=500)
        smoothed = smooth(raw_values, window='hanning', window_len=200)
        self.assertTrue(isinstance(smoothed, np.ndarray))
        self.assertEqual(len(smoothed), 500)
        self.assertEqual(round(smoothed[5], 2), round(smoothed[10], 2))
        self.assertNotEqual(round(raw_values[5], ), round(raw_values[10], 2))

