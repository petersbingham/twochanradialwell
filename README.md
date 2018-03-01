# TwoChanRadialWell
Calculates solutions to the two channel radial well as described in Newton's "Scattering Theory of Waves and Particles".

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/TwoChanRadialWell.git
    cd TwoChanRadialWell
    python setup.py install
    
## Dependencies
Standard Libraries: numpy & scipy

pynumwrap https://github.com/petersbingham/pynumwrap

pynumutil https://github.com/petersbingham/pynumutil

channelutil https://github.com/petersbingham/channelutil

## Usage

The getSmatFun function returns a function reference to the S-matrix as a function of energy. It's signiture looks like:
```python
getSmatFun(r0, v1, v2, chanCalc, lam)
```
`r0`, `v1`, `v2` and `lam` should be obvious after consulting Newton's text. The channel calc is created in the client code and is described at the link in the Dependencies section. It contains the threshold values.

There are two types that TwoChanRadialWell is compatible with, standard python types and mpmath types. Python types is the default. To change to mpmath types call the module function `usempmathTypes()`. The example below illustrate usage.
```python
>>> import TwoChanRadialWell as radwell
>>> import channelutil as chanutil
>>> chanCalc = chanutil.calculator([0.,2.], massMult=chanutil.MASSMULT_HARTREES)
>>> sMatfFn = radwell.getSmatFun(1., 2., 2., chanCalc, 1.)
>>> print sMatfFn(0.)
[[  1.00000000+0.j   0.00000000+0.j]
 [  0.00000000+0.j -13.56891277+0.j]]
```
