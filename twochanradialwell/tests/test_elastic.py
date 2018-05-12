import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import twochanradialwell as radwell
import channelutil as chanutil
import unittest

class test_elastic(unittest.TestCase):
    def runTest(self):
        asymcalc = chanutil.AsymCalc(chanutil.HARTs, [0,0])
        fun = radwell.getSmatFun(1., 2., 2., asymcalc, 1.)

        expectMat = radwell.nw.matrix([[1., 0.],[0., 1.]])
        gotMat = fun(0.)
        self.assertTrue(radwell.nw.np.allclose(gotMat,expectMat))

        expectMat = radwell.nw.matrix(\
          [[-0.31507906+0.89937359j, -0.28602697-0.10020431j],
           [-0.28602697-0.10020431j, -0.31507906+0.89937359j]])
        gotMat = fun(1.8)
        self.assertTrue(radwell.nw.np.allclose(gotMat,expectMat))

        expectMat = radwell.nw.matrix(\
          [[0.1665721+0.95554805j, -0.23965872+0.04177755j],
           [-0.23965872+0.04177755j, 0.1665721+0.95554805j]])
        gotMat = fun(3.0)
        self.assertTrue(radwell.nw.np.allclose(gotMat,expectMat))



if __name__ == "__main__":
    #Just for debug
    b = test_elastic()
    b.runTest()
