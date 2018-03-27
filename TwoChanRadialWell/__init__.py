import numpy as np
import scipy.linalg as la

import pynumutil as num
import pynumwrap as nw
try:
    import tisutil as tu
except:
    tu = None

EQUIVALENT_TESTS = False
LIN_ALGEBRA = False

RESULTS_TYPE_FLOAT32 = -1
RESULTS_TYPE_FLOAT64 = -2
RESULTS_TYPE_DEFAULT = 0

gu = num.nearlyEqual()

class RadWellException(Exception):
    def __init__(self, string):
        self.string = string
    def __str__(self):
        return "Rad Well Error: " + self.string

class mat:
    PADDING = 3
    def __init__(self, size, precision):
        self.size = size
        self.precision = precision
        self.min = pow(10,-self.precision)

    def __getitem__(self, i):
        return self._getRow(i)

    def getMatrix(self):
        mlist = []
        for m in range(self.size):
            rlist = []
            for n in range(self.size):
                rlist.append(nw.complex(self[m][n]))
            mlist.append(rlist)
        return nw.matrix(mlist)

    def __str__(self):
        isImag = self._isImag()

        maxLen = 0
        for m in range(self.size-1):
            for n in range(self.size-1):
                eLen = len(self._getFormattedStr(self[m][n], isImag))
                if eLen > maxLen:
                    maxLen = eLen
        s = ""
        for m in range(self.size):
            for n in range(self.size):
                s += self._padStr(self[m][n], isImag, maxLen)
            s += "\n"
        return s

    def _padStr(self, value, isImag, maxLen=None):
        theStr = self._getFormattedStr(value, isImag)
        if maxLen is not None:
            return theStr.ljust(maxLen + self.PADDING)
        else:
            return theStr

    def _isImag(self):
        for m in range(0,self.size):
            for n in range(0,self.size):
                if abs(float(nw.complex(self[m][n]).imag)) > self.min:
                    return True
        return False

    def _getFormattedStr(self, value, isImag):
        if isImag:
            return nw.formattedComplexString(nw.complex(value), self.precision)
        else:
            return nw.formattedFloatString(nw.complex(value).real, self.precision)

class Mats:
    def __init__(self, v1, v2, asymCal, lam):
        self.K = Mats.Kmat(asymCal)
        self.V = Mats.Vmat(v1, v2, asymCal, lam)
        self.A = Mats.Amat(self.K, self.V)
        if LIN_ALGEBRA:
            self.aSq = Mats.aSqMat(self.A)
            self.a = Mats.aMat(self.aSq)

    def setEnergy(self, ene):
        self.K.setEnergy(ene)
        if LIN_ALGEBRA:
            self.aSq.calculate()
            self.a.calculate()

    def printMats(self):
        print "\nK:"
        print str(self.K) + "\n\nV:"
        print str(self.V) + "\n\nA:"
        print str(self.A) + "\n\nsqrt(A):"
        print str(self.aSq) + "\n\na:"
        print str(self.a)

    class Kmat(mat):
        def __init__(self, asymCal):
            mat.__init__(self, 2, num.precision)
            self.ene = 0
            self.asymCal = asymCal
        def setEnergy(self, ene):
            self.ene = ene
        def _getRow(self, i):
            if i==0:
                return [self.k(0), 0]
            else:
                return [0, self.k(1)]
        def k(self, ch):
            return self.asymCal.k(ch, self.ene)

    class Vmat(mat):
        def __init__(self, v1, v2, asymCal, lam):
            mat.__init__(self, 2, num.precision)
            self.v1 = v1
            self.v2 = v2
            self.massMult = asymCal.getEneConv()
            self.lam = lam
        def _getRow(self, i):
            if i==0:
                return [-self.v1*self.massMult, -0.5*self.lam*self.massMult]
            else:
                return [-0.5*self.lam*self.massMult, -self.v2*self.massMult]

    class Amat(mat):
        def __init__(self, K, V):
            mat.__init__(self, 2, num.precision)
            self.K = K
            self.V = V
        def _getRow(self, i):
            return [nw.pow(self.K[i][0],2)-self.V[i][0],
                    nw.pow(self.K[i][1],2)-self.V[i][1]]

    class aSqMat(mat):
        def __init__(self, A):
            mat.__init__(self, 2, num.precision)
            self.A = A
            self.calculate()
        def calculate(self):
            m = self.A.getMatrix()
            eigvals,eigvecs = np.linalg.eig(m)
            self.aSq = np.diagflat(eigvals) 
        def _getRow(self, i):
            return [self.aSq[i,0], self.aSq[i,1]]

    class aMat(mat):
        def __init__(self, aSq):
            mat.__init__(self, 2, num.precision)
            self.aSq = aSq
            self.calculate()
        def calculate(self):
            m = self.aSq.getMatrix()
            self.a = la.sqrtm(m)
        def _getRow(self, i):
            return [self.a[i,0], self.a[i,1]]

class Smat(mat):
    def __init__(self, r0, mats, resultsType=RESULTS_TYPE_DEFAULT):
        mat.__init__(self, 2, num.precision)
        self.resultsType = resultsType
        self.numChannels = 2
        self.mats = mats
        self.r0 = r0

    def setEnergy(self, ene):
        self.mats.setEnergy(ene)
        return self

    def _getRow(self, i):
        if i==0:
            return [self._S_11(), self._S_12()]
        else:
            return [self._S_21(), self._S_22()]

#######

    def _truncate_float(self, num, digits):
        s = "{:."+str(digits)+"e}"
        m, e = s.format(num).split('e')
        return float(float(m)*float(10**float(e)))

    def _cast_result(self, result):
        if self.resultsType != RESULTS_TYPE_DEFAULT:
            if self.resultsType == RESULTS_TYPE_FLOAT64:
                return complex(result)
            elif self.resultsType == RESULTS_TYPE_FLOAT32:
                return complex(float(np.float32(result.real)), 
                               float(np.float32(result.imag)))
            else:
                return complex(self._truncate_float(float(result.real), 
                                                    self.resultsType), 
                               self._truncate_float(float(result.imag), 
                                                    self.resultsType))
        else:
            return result

    def _S_11(self):
        ret = self._g(-self._rho_1(), self._rho_2()) / self._denum() *\
                self._exp(2.0*self._rho_1())
        return self._cast_result(ret)

    def _S_12(self):
        ret = 2.0 * (self._zeta_1()-self._zeta_2()) *\
            nw.sqrt(self._alp_1()*self._alp_2()*self._rho_1()*self._rho_2()) /\
            self._denum() *\
            self._exp(self._rho_1() + self._rho_2())
        return self._cast_result(ret)

    def _S_21(self):
        ret = self._S_12()
        return self._cast_result(ret)

    def _S_22(self):
        ret = self._g(self._rho_1(), -self._rho_2()) / self._denum() *\
        self._exp(2.0*self._rho_2())
        return self._cast_result(ret)

    def _g(self, rho_1, rho_2):
        complex1 = rho_1 *\
                   (self._zeta_1()*self._alp_1() - self._zeta_2()*self._alp_2())
        complex2 = rho_2 *\
                   (self._zeta_2()*self._alp_1() - self._zeta_1()*self._alp_2())
        real = (self._alp_1()-self._alp_2()) *\
               (rho_1*rho_2 - self._zeta_1()*self._zeta_2())
        return real + (complex1+complex2)*1.0j

    def _denum(self, test=True):
        value = self._g(self._rho_1(), self._rho_2())
        if test and gu.complexCompare(value, 0.0):
            raise RadWellException("_denum: Zero")
        return value

    def _exp(self, rho):
        return nw.exp(-1.0j*rho)

#######

    def _rho_1(self):
        return self._rho_alp(0)

    def _rho_2(self):
        return self._rho_alp(1)

    def _alp_1(self):
        cal1 = nw.pow(self._e_n(0), 2) - nw.pow(self._R_alp(0), 2)
        cal2 = nw.pow(self._R_alp(1), 2) - nw.pow(self._e_n(1), 2)
        if EQUIVALENT_TESTS:
            if not gu.complexCompare(cal1, cal2):
                raise RadWellException("_alp_1: " + str(cal1) + "   " + str(cal2))
        return cal1

    def _alp_2(self):
        cal1 = nw.pow(self._e_n(1), 2) - nw.pow(self._R_alp(0), 2)
        cal2 = nw.pow(self._R_alp(1), 2) - nw.pow(self._e_n(0), 2)
        if EQUIVALENT_TESTS:
            if not gu.complexCompare(cal1, cal2):
                raise RadWellException("_alp_2: " + str(cal1) + "   " + str(cal2))
        return cal1

    def _zeta_1(self):
        return self._e_n(0) / nw.tan(self._e_n(0))

    def _zeta_2(self):
        return self._e_n(1) / nw.tan(self._e_n(1))

#######

    def _e_n(self, n):
        if LIN_ALGEBRA:
            cal1 = self.mats.a[n][n] * self.r0
        cal2 = nw.sqrt(self._e_n_Sq_alt(n))
        if EQUIVALENT_TESTS:
            if not gu.complexCompare(cal1, cal2):
                raise RadWellException("_e_n: " + str(cal1) + "   " + str(cal2))
        if LIN_ALGEBRA:
            return cal1
        else:
            return cal2

    def _e_n_Sq_alt(self, n):
        first = (nw.pow(self._R_alp(0),2.0)+nw.pow(self._R_alp(1),2.0)) / 2.0
        a = nw.pow(nw.pow(self._R_alp(0),2.0)-nw.pow(self._R_alp(1),2.0),2.0)
        b = 4.0*nw.pow(self.mats.V[0][1],2.0)*nw.pow(self.r0,4.0)
        second = nw.sqrt( a + b ) / 2.0
        if n==0:
            return first+second
        else:
            return first-second

    def _rho_alp(self, ch):
        return self.mats.K.k(ch) * self.r0

    def _R_alp(self, ch):
        cal1 = nw.sqrt( self.mats.A[ch][ch] ) * self.r0
        cal2 = nw.sqrt( nw.pow(self._rho_alp(ch),2) -\
        self.mats.V[ch][ch]*nw.pow(self.r0,2.0) )
        if EQUIVALENT_TESTS:
            if not gu.complexCompare(cal1, cal2):
                raise RadWellException("_R_alp: " + str(cal1) + "   " + str(cal2))
        return cal2

########################################################################   
######################### Public Interface #############################
########################################################################

def getSmatFun(r0, v1, v2, asymCal, lam, resultsType=RESULTS_TYPE_DEFAULT):
    mats = Mats(v1, v2, asymCal, lam)
    sMat = Smat(r0, mats, resultsType)
    funPtr = lambda ene : sMat.setEnergy(ene).getMatrix()
    if tu is not None:
        assert asymCal.getUnits() == tu.HARTs
        return tu.cSmat(funPtr, asymCal)
    else:
        return funPtr

def usePythonTypes(dps=nw.dps_default_python):
    nw.usePythonTypes(dps)

def useMpmathTypes(dps=nw.dps_default_mpmath):
    nw.useMpmathTypes(dps)

def setTypeMode(mode, dps=None):
    nw.setTypeMode(mode, dps)
