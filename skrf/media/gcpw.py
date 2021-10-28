

"""
gcpw (:mod:`skrf.media.gcpw`)
========================================

.. autosummary::
   :toctree: generated/

   GCPW

"""
from scipy.constants import epsilon_0, mu_0, c
from scipy.special import ellipk
from numpy import pi, sqrt, log, tanh, exp, zeros, ones
from .media import Media
from ..tlineFunctions import surface_resistivity
from ..constants import NumberLike
from typing import Union, TYPE_CHECKING

if TYPE_CHECKING:
    from .. frequency import Frequency


class GCPW(Media):
    """
    Coplanar waveguide over ground

    This class was derived from the technical documentation [#]_ provided
    by the qucs project [#]_ .
    The variables and properties of this class are coincident with
    their derivations.

    Parameters
    ----------
    frequency : :class:`~skrf.frequency.Frequency` object, optional
        frequency band of the media. The default is None.
    z0 : number or array-like, optional
        the port impedance for media. The default is None.
        Only needed if it's different from the characteristic impedance
        of the transmission.
    w : number or array-like, optional
        width of center conductor, in m. Default is 2e-3.
    s : number or array-like
        width of gap, in m. Default is 5e-4.
    h : number or array-like
        height of substrate, in m. Default is 1.5e-3.
    ep_r : number or array-like, optional
        relative permittivity of substrate. Default is 4.5.
    t : number, or array-like, optional
        conductor thickness, in m. Default is None (metalization thickness neglected)
    rho : number, or array-like, optional
        resistivity of conductor. Default is None (conductor losses neglected)

    References
    ----------
    .. [#] http://qucs.sourceforge.net/docs/technical.pdf
    .. [#] http://www.qucs.sourceforge.net/

    """
    def __init__(self, 
                frequency: Union['Frequency', None] = None,
                z0: Union[NumberLike, None] = None,
                w: NumberLike = .002, 
                s: NumberLike = .0005,
                h: NumberLike = .0015,
                ep_r: NumberLike = 4.5,
                t: Union[NumberLike, None] = None,
                rho: Union[NumberLike, None] = None,
                *args, **kwargs):

        Media.__init__(self, frequency=frequency,z0=z0)

        self.w, self.s, self.h, self.ep_r, self.t, self.rho =\
                w, s, h, ep_r, t, rho


    def __str__(self) -> str:
        f=self.frequency
        output =  \
                'Grouned Coplanar Waveguide Media.  %i-%i %s.  %i points'%\
                (f.f_scaled[0],f.f_scaled[-1],f.unit, f.npoints) + \
                '\n W= %.2em, S= %.2em'% \
                (self.w,self.s)
        return output


    def __repr__(self) -> str:
        return self.__str__()


    @property
    def ep_re0(self) -> NumberLike:
        """
        Effective permittivity of the grounded coplanar waveguide at f=0.

        See equations 12.10, 12.11, and 12.18 in the qucs documentation.
        """
        kr1 = self.K_ratio(self.k1)
        kr3 = self.K_ratio(self.k3)
        q = kr3 / (kr1 + kr3)

        er = 1 + q * (self.ep_r - 1)

        if self.t is None:
            return er

        # Compensate for metallization thickness
        tr = 0.7 * self.t / self.s

        return er - (tr * (er - 1) / (kr1 + tr))


    @property
    def ep_re(self) -> NumberLike:
        """
        Effective permittivity of the grounded coplanar waveguide.
        """
        w, s, h, ep_r = self.w, self.s, self.h, self.ep_r

        p = log(w / h)
        u = 0.54 - 0.64 * p + 0.015 * p * p
        v = 0.43 - 0.86 * p + 0.54 * p * p
        g = exp(u * log(w / s) + v)

        fte = c / (4 * h * sqrt(ep_r - 1))
        ser0 = sqrt(self.ep_re0)

        num = sqrt(ep_r) - ser0
        den = 1 + g * (self.frequency.f / fte) ** -1.8

        return (ser0 + (num / den)) ** 2


    @property
    def k1(self) -> NumberLike:
        """
        Intermediary parameter. see qucs docs on cpw lines.

        Defined as:

        .. math::

                k = \\frac{w}{w + 2s}

        """
        w, s = self.w, self.s
        return w / (w + s + s)


    @property
    def ke(self) -> NumberLike:
        """
        Intermediary parameter used in correcting for metalization thickness.

        See equation 12.17 in the qucs documentation.
        """
        t, w, s, k1 = self.t, self.w, self.s, self.k1

        d = 1.25 * t * (1 + log(4 * pi * w / t)) / pi
        return k1 + d * (1 - (k1 ** 2)) / (2 * s)


    @property
    def k3(self) -> NumberLike:
        """
        Intermediary parameter.

        See equation 12.12 in the qucs documentation.

        """
        w, h, s = self.w, self.h, self.s

        f = lambda x: tanh((pi * x) / (h * 4))

        return f(w) / f(w + 2 * s)


    @property
    def Z0(self) -> NumberLike:
        """
        Characteristic impedance of the transmission line.

        See equation 12.13 in the qucs documentation.
        """
        
        # Compensate for metallization thickness
        if self.t is None:
            kr1 = self.K_ratio(self.k1)
        else:
            kr1 = self.K_ratio(self.ke)

        kr3 = self.K_ratio(self.k3)

        return (60 * pi) / (sqrt(self.ep_re) * (kr1 + kr3))


    @property
    def alpha_conductor(self) -> NumberLike:
        """
        Losses due to conductor resistivity.

        See equation 12.26 in the qucs documentation.

        Returns
        -------
        alpha_conductor : array-like
                attenuation due to conductor losses

        See Also
        --------
        surface_resistivity : calculates surface resistivity
        """
        if self.rho is None:
            raise AttributeError('Conductor loss calculation requires conductivity to be defined.')

        if self.t is None:
            raise AttributeError('Conductor loss calculation requires metallization thickness to be defined.')

        t, s, w, k1, ep_re = self.t, self.s, self.w, self.k1, self.ep_re
        r_s = surface_resistivity(f=self.frequency.f, rho=self.rho, mu_r=1)

        a = w / 2
        b = a + s
        ks = 1 - (k1 ** 2)

        K = ellipk(k1)
        K_p = ellipk(sqrt(ks))

        term = lambda x: (pi + log(8 * pi * x * (1 - k1) / (t * (1 + k1)))) / x

        return (r_s * sqrt(ep_re) * (term(a) + term(b))) / (480 * pi * K * K_p * ks)


    @property
    def alpha_dielectric(self) -> NumberLike:
        """
        Losses in dielectric.

        See equation 11.79 in the qucs documentation.

        Returns
        -------
        alpha_dielectric : array-like
                attenuation due to conductor losses
        """
        pass
    
 
    @property
    def gamma(self) -> NumberLike:
        """
        Propagation constant.

        See Also
        --------
        alpha_conductor : calculates attenuation due to conductor losses
        """
        beta = 1j * 2 * pi * self.frequency.f * sqrt(self.ep_re * epsilon_0 * mu_0)

        if (self.rho is None) or (self.t is None):
            return beta

        return beta + self.alpha_conductor


    @staticmethod
    def K_ratio(k: NumberLike) -> NumberLike:
        """
        K_ratio is the approximate ratio of the complete elliptic integral of k to its complement.

        See equations 12.4 and 12.5 in qucs documentation.
        """

        if k < 0 or k > 1:
            raise ValueError("Function defined on range 0 <= k <= 1")

        term = lambda x: log(2 * (1 + sqrt(x)) / (1 - sqrt(x)))

        if k > (1 / sqrt(2.0)):
            return term(k) / pi

        k_p = sqrt(1 - k ** 2)
        return pi / term(k_p)

