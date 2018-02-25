import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import TwoChanRadialWell as radwell

import unittest

class test_elastic(unittest.TestCase):
    def runTest(self):
        self.assertEqual(True,True)

if __name__ == "__main__":
    #Just for debug
    b = test_elastic()
    b.runTest()
