'''
.. module:: skrf.media.gcpw

========================================
cpw (:mod:`skrf.media.gcpw`)
========================================

Coplanar waveguide over ground class

This class was made from the technical documentation [#]_ provided
by the qucs project [#]_ .
The variables  and properties of this class are coincident with
their derivations.

.. [#] http://qucs.sourceforge.net/docs/technical.pdf
.. [#] http://www.qucs.sourceforge.net/

.. autosummary::
   :toctree: generated/

   GCPW

'''

from numpy import pi,sqrt,log,ones,tanh,piecewise,nan,exp
from scipy.special import ellipk
from scipy.constants import c, epsilon_0, mu_0
from .cpw import CPW

class GCPW(CPW):
    '''
    Coplanar Waveguide Over Ground Media

    Parameters
    -------------
    frequency : :class:`~skrf.frequency.Frequency` object, optional
        frequency band of the media. The default is None.
    z0 : number, array-like, optional
        the port impedance for media. The default is None.
        Only needed if it's different from the characterisitc impedance 
        of the transmission line.
    w : number, or array-like, optional
            width of center conductor, in m. Default is 70.
    s : number, or array-like, optional
            width of gap, in m. Default is 4.
    h : number, or array-like, optional
            height of substrate, in m. Default is 1.
    ep_r : number, or array-like, optional
            relative permittivity of substrate. Default is 3.
    tan_delta : number, or array-like, optional
        loss tangent of substrate. Default is 0 (lossless medium) 
    t : number, or array-like, optional
            conductor thickness, in m. Default is None (metallization thickness neglected)
    rho: number, or array-like, optional
            resistivity of conductor. Default is None (resistivity neglected)
    '''
    def __init__(self, frequency=None, z0=None, w=70, s=4, h=1,
                 ep_r=3, tan_delta=0, t=None, rho=None,  *args, **kwargs):

        super().__init__(frequency=frequency, z0=z0, w=w, s=s,
            ep_r=ep_r, t=t, rho=rho)

        self.h = h    
        self.tan_delta = tan_delta

    @property
    def ep_re0(self): # tested - correct
        '''
        Effective permittivity of the grounded CPW without dispersion.

        Calculated per equations 12.10 and 12.18 of Qucs technical paper.
        '''
        e0 = 1 + self.q * (self.ep_r - 1)

        if self.t is None:
            return e0

        ts = 0.7 * self.t / self.s

        num = (e0 - 1) * ts
        den = _K_over_Kprime(self.k1) + ts

        return e0 - (num / den)
    
    @property
    def q(self):
        '''
        Filling factor per equation 12.11 of Qucs technical paper.
        '''
        q3 = _K_over_Kprime(self.k3)
        return q3 / (_K_over_Kprime(self.k1) + q3)

    @property
    def k3(self):
        '''
        Intermediary parameter - per equation 12.12 of Qucs technical paper.
        '''
        t = lambda x: tanh((pi * x) / (4 * self.h))
        return t(self.w) / t(self.w + (2 * self.s))

    @property
    def ke(self):
        '''
        Intermediary parameter - per equation 12.17 of Qucs technical paper.
        '''
        delta = (1.25 * self.t / pi) * (1 + log((4 * pi * self.w) / self.t))
        k1 = self.k1

        return k1 + (1 - (k1 ** 2)) * delta / (2 * self.s)

    @property
    def ep_re(self): # tested - correct
        '''
        Effective permittivity of the grounded CPW.

        Calculated per equation 12.19 of Qucs technical paper.
        '''
        sr_er0 = sqrt(self.ep_re0)

        num = sqrt(self.ep_r) - sr_er0
        den = 1 + (self. G * (self.frequency.f / self.f_te) ** (-1.8))

        return (sr_er0 + (num / den)) ** 2

    @property
    def G(self): # tested - correct
        '''
        Intermediary parameter - per equation 12.20 of Qucs technical paper.
        '''
        p = log(self.w / self.h) # equation 12.23
        u = 0.54 - 0.64 * p + 0.015 * (p ** 2) # equation 12.21
        v = 0.43 - 0.86 * p + 0.54 * (p ** 2) # equation 12.22

        return exp(u * log(self.w / self.s) + v)

    @property
    def f_te(self):
        '''
        TE0 cutoff frequency per equation 12.24 of Qucs technical paper.
        '''
        return c / (4 * self.h * sqrt(self.ep_r - 1)) 

    @property
    def Z0(self): # close to correct when using approximation for K/K'
        '''
        Characteristic impedance.

        Calculated per equation 12.13 of Qucs technical paper.
        '''

        if self.t is None:
            k1 = self.k1
        else:
            k1 = self.ke

        z = 60 * pi
        z = z / (_K_over_Kprime(k1) + _K_over_Kprime(self.k3))
        return z / sqrt(self.ep_re)

    @property
    def alpha_dielectric(self):
        '''
        Dielectric loss factor per equation 11.79 of Qucs technical paper.
        '''
        ep_re = self.ep_re
        return (self.ep_r / (self.ep_r - 1)) * ((ep_re - 1) / sqrt(ep_re)) * pi * self.tan_delta / c

    @property
    def gamma(self):
        '''
        Propagation constant


        See Also
        --------
        alpha_conductor : calculates conductor losses
        alpha_dielectric : calculates dielectric losses
        '''
        alpha = self.alpha_dielectric

        if self.rho is not None and self.t is not None:
            alpha += self.alpha_conductor

        #beta = 1j * 2 * pi * self.frequency.f * sqrt(self.ep_re * epsilon_0 * mu_0)

        beta = 1j * self.frequency.f * sqrt(self.ep_re) * 2 * pi / c

        return alpha + beta

# Supporting functions
def _K_over_Kprime_approx(k):

    conditions = [
        k >= 0 and k < (1/sqrt(2)),
        k >= (1/sqrt(2)) and k < 1,
    ]

    functions = [
        lambda x: pi / log(2 * (1 + sqrt(_comp_mod(x))) / (1 - sqrt(_comp_mod(x)))),
        lambda x: log(2 * (1 + sqrt(x)) / (1 - sqrt(x))) / pi,
        nan,
    ]
    
    return piecewise(k, conditions, functions)

def _K_over_Kprime(k):
    return _K_over_Kprime_approx(k)
    #return ellipk(k) / ellipk(_comp_mod(k))

def _comp_mod(k):
    return sqrt(1 - k**2)