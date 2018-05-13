# twochanradialwell
Calculates solutions to the two channel radial well as described in Newton's "Scattering Theory of Waves and Particles".

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/twochanradialwell.git
    cd twochanradialwell
    python setup.py install
    
## Dependencies
Standard Libraries: 
 - numpy
 - scipy

Author Libraries (these will have their own dependencies):
 - pynumwrap https://github.com/petersbingham/pynumwrap
 - pynumutil https://github.com/petersbingham/pynumutil
 - channelutil https://github.com/petersbingham/channelutil

## Usage

The getSmatFun function returns a function reference to the S-matrix as a function of energy. It's signature looks like:
```python
get_Smat_fun(r0, v1, v2, asymcalc, lam)
```
`r0`, `v1`, `v2` and `lam` should be obvious after consulting Newton's text. The channel calc is created in the client code and is described at the link in the Dependencies section. It contains the threshold values.

The example below illustrates usage.
```python
>>> import twochanradialwell as radwell
>>> import channelutil as chanutil
>>> asymcalc = chanutil.AsymCalc(chanutil.hartrees, thresholds=[0.,2.])
>>> smatfun = radwell.get_Smat_fun(1., 2., 2., asymcalc, 1.)
>>> print smatfun(0.)
[[  1.00000000+0.j   0.00000000+0.j]
 [  0.00000000+0.j -13.56891277+0.j]]
```

There are two types that twochanradialwell is compatible with, standard python types and mpmath types. Python types is the default. To change to mpmath types call the module function `use_mpmath_types()`.
