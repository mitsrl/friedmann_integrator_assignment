"""Friedmann equation integrators"""


import numpy as np
from scipy.integrate import cumtrapz


# We are going to use a *class* for our Friedmann integrator. Classes are used
# for many things, but here I'm using one to closely link *functionality*
# with *data*. In functional programming you have functions that operate
# on data (variables/parameters) to produce new data. In so-called object
# oriented programming, classes/objects internally contain both data and functionality.
# This is logically useful, since most functions only work given valid input
# data, and classes provide a way to ensure validity and proper use.


# Our integrator class contains two types of data: the input cosmological
# parameters, and the results of the numerical integration. It contains functions
# (a function that is part of a class is called a "method"). Some of these
# methods produce one type of data (the integration) from another (the
# cosmological parameters). Other methods manipulate data into a different
# form, converting comoving distance to luminocity distance (for example).

class CumulativeSumIntegrator:
    """Friedmann equation integrator using a cumulative Riemann sum.

    Parameters
    ----------
    Omega_r: float between 0 and 1
        Radiation density as a fration of the critical density.
    Omega_m: float between 0 and 1
        Matter density as a fration of the critical density.
        Omega_r + Omega_m should be less than or equal to 1.
    H0: float > 0
        Hubble constant in units of km/s/Mpc.

    """

    def __init__(self, Omega_r, Omega_m, H0):
        # *self* is an implicit method parameter. It is how you refer to the
        # object from within the object.
        Omega_lam = 1 - Omega_r - Omega_m

        # Always good to check user input.
        if H0 <= 0:
            raise ValueError("Hubble constant must be positive.")
        for om in (Omega_m, Omega_r, Omega_lam):
            if om < 0 or om > 1:
                raise ValueError("Invalid Omegas.")

        # Here we *store* data in the object. We are renaming the variables to begin
        # with an '_'. This is called *encapsulation*, and signals to users of the
        # code not to change thier values or risk getting inconsistent results.
        self._Omega_r = Omega_r
        self._Omega_m = Omega_m
        self._Omega_lam = Omega_lam
        self._H0 = H0

        # These are wrong! Should be calculated from H0 and 
        # put them in to the correct units of years and Mpc.
        self._hubble_time = 1 / H0
        self._hubble_distance = 2997.9246

        # Choose scale factors at which we will evaluate the integrals. I'm choosing
        # logarithmic spacing here, since the P-set askes for a log-log plot.
        # The last point in this array is a=1, the value at the present day.
        self._a = np.logspace(-7, 0, 1000, endpoint=True)

        # Most of the setup for the class is complete. Perform the numerical integrals.
        self._integrate()

    def _integrate(self):
        """This is wrong unless Omega_m = 1."""

        # Time variables multiplied by H0, which is unitless.
        # In general in numerics you want to be dealing with
        # dimensionless numbers, adding in dimensionful factors after the calculation.
        self._tH0 = 2 / 3 *self._a ** (3 / 2)
        self._etaH0 = 2 * np.sqrt(self._a)

        # To implement this propertly, write down the differential equations for eta and t
        # in terms of the scale factor a. Notice that these are separable so the solution may
        # be written as an integral. Unfortunately, the integral is not analytic so must
        # be done numerically.  Write an approximation to the integrals as a Riemann sum (with
        # non-uniform step spacing). Code up the sums, making use of the functions numpy.diff() and
        # numpy.cumsum() (see the online documentation for these functions).

    # These *properties* provide users of this class *read-only* access to encapsulated data and
    # simple transformations thereof.
    @property
    def Omega_r(self):
        return self._Omega_r

    @property
    def Omega_m(self):
        return self._Omega_m

    @property
    def Omega_lam(self):
        return self._Omega_lam

    @property
    def H0(self):
        return self._H0

    @property
    def a(self):
        return self._a.copy()

    @property
    def z(self):
        return 1 / self._a - 1

    @property
    def hubble_parameter(self):
        """This is wrong unless Omega_m = 1."""
        return self.H0 / self._a**2

    @property
    def comoving_distance(self):
        return self._hubble_distance * (self._etaH0[-1] - self._etaH0)

    @property
    def luminocity_distance(self):
        raise NotImplementedError("Your job.")


# One of the advantages of object oriented programming is extensibility. You can
# make a subclass, which is just a copy of a class the changes (overwrites) just
# a few of the methods. Reimplement the `integrate` method in the class below
# using scipy.integrate.cumtrapz. If you end up with code that is repeated in both
# implementations, you should move them to a different method so you don't have to
# repeated code. Compare the accuracy of the two integrators.

class TrapzIntegrator(CumulativeSumIntegrator):

    def integrate(self):
        raise NotImplementedError("Your Job")
