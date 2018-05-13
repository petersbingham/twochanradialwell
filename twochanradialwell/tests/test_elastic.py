import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import twochanradialwell as radwell
import channelutil as chanutil
import unittest

class test_elastic(unittest.TestCase):
    def runTest(self):
        asymcalc = chanutil.AsymCalc(chanutil.hartrees, [0,0])
        fun = radwell.get_Smat_fun(1., 2., 2., asymcalc, 1.)

        expect_mat = radwell.nw.matrix([[1., 0.],[0., 1.]])
        got_mat = fun(0.)
        self.assertTrue(radwell.nw.np.allclose(got_mat,expect_mat))

        expect_mat = radwell.nw.matrix(\
          [[-0.31507906+0.89937359j, -0.28602697-0.10020431j],
           [-0.28602697-0.10020431j, -0.31507906+0.89937359j]])
        got_mat = fun(1.8)
        self.assertTrue(radwell.nw.np.allclose(got_mat,expect_mat))

        expect_mat = radwell.nw.matrix(\
          [[0.1665721+0.95554805j, -0.23965872+0.04177755j],
           [-0.23965872+0.04177755j, 0.1665721+0.95554805j]])
        got_mat = fun(3.0)
        self.assertTrue(radwell.nw.np.allclose(got_mat,expect_mat))



if __name__ == "__main__":
    #Just for debug
    b = test_elastic()
    b.runTest()
