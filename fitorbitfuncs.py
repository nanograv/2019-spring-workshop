#! /usr/bin/python

from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import sys

d2r = np.pi / 180

def mean_anomaly(pb, t, t0, pbdot=0):
    """
    Computes mean anomaly, given orbital period and time parameters.

    Inputs:
        - pb = orbital period [days]
        - t = epoch to evaluate angle [MJD]
        - t0 = reference epoch [MJD]
        - pbdot = time derivative of orbital period [  ]

    Output:
        - mean anomaly [deg]
    """

    # check if input time is a list or NumPy array.
    if (isinstance(t, list)):
        t = np.array(t)
    elif (isinstance(t, np.ndarray)):
        pass

    # now do calculation
    dt = t - t0
    pbdot *= 86400.
    ma = 360 / pb * (dt - 0.5 * pbdot / pb * dt**2) % 360

    # make sure that 0 < ma < 360
    if (isinstance(ma, np.ndarray) and np.any(ma < 0)):
        ma[np.where(ma < 0)] += 360
    elif (not isinstance(ma, np.ndarray) and ma < 0):
        ma += 360

    return ma

def ecc_anomaly(ma, ecc, ea0=0.5):
    """
    Computes eccentric anomaly, given mean anomaly and eccentricity.

    Inputs:
        - ma = mean anomaly [deg]
        - ecc = orbital eccentricity [  ]
        - ea0 = initial guess of eccentric anomaly [  ]

    Output:
        - eccentric anomaly [deg]
    """

    ma_in = ma * d2r
    ea = 0

    # if MA is an array, loop over each entry and apply N-R method.
    if (isinstance(ma, np.ndarray)):
       count = 0
       ea = np.zeros(len(ma))

       # in this case, turn 'ecc' into an array.
       if (not isinstance(ecc, np.ndarray)):
           ecc = np.zeros(len(ma)) + ecc

       # compute EA for each MA, separately.
       for ma0, ecc0 in zip(ma_in, ecc):
           ea_mid = ma0
           for i in range(100):
               f  = ea_mid - ecc0 * np.sin(ea_mid) - ma0
               fp = 1 - ecc0 * np.cos(ea_mid)
               ea_mid -= f / fp
               if (np.fabs(ea_mid - ea0) < 1e-12):
                   ea[count] = ea_mid
                   count += 1
                   break
               ea0 = ea_mid

    # otherwise, do single calculation and leave as scalar.
    else:
        ea = ma
        for i in range(100):
           f  = ea - ecc * np.sin(ea) - ma
           fp = 1 - ecc * np.cos(ea)
           ea -= f / fp
           if (np.fabs(ea - ea0) < 1e-12):
               break
           ea0 = ea

    ea /= d2r

    # make sure that 0 < EA < 360
    if (isinstance(ea, np.ndarray) and np.any(ea < 0)):
        ea[np.where(ea < 0)] += 360
    elif (not isinstance(ea, np.ndarray) and ea < 0):
        ea += 360

    return ea

def true_anomaly(ea, ecc):
    """
    Computes true anomaly, given eccentric anomaly and eccentricity.

    Inputs:
        - ea = eccentric anomaly [deg]
        - ecc = orbital eccentricity [deg]

    Output:
        - true anomaly [deg]
    """

    ea_in = ea * d2r
    ta = 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ea_in / 2)) / d2r

    # make sure that 0 < TA < 360
    if (isinstance(ta, np.ndarray) and np.any(ta < 0)):
        ta[np.where(ta < 0)] += 360
    elif (not isinstance(ta, np.ndarray) and ta < 0):
        ta += 360

    return ta

def peri_argument(om0, pb, ta, omdot=0):
    """
    Computes periastron argument as a function of time, given initial 
    value, orbital pb, true anomaly (ta) and periastron advance.

    Inputs:
        - om0 = periastron argument [deg]
        - pb = orbital period [days]
        - ta = true anomaly [deg]
        - omdot = time derivative of periastron argument [deg / yr]

    Output:
        - periastron argument [deg]
    """
    pb_in = pb / 365.25 # convert to years.
    om = om0 + omdot * ta * pb / 2 / np.pi % 360
    return om

def mass_function(x, pb):
    """
    Computes Keplerian mass function, given projected size and orbital period.

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]

    Output:
        - mass function [solar mass]
    """

    pb_in = pb * 86400
    T_sun = 4.925490947 * 1e-6

    return 4 * np.pi**2 * x**3 / T_sun / pb_in**2

def m2_massfunction(mp, sini, mf, n_iter=100):
    """
    Computes the value of the companion mass based on the mass function, 
    for given/fixed values of the pulsar mass and sine of inclination.

    Inputs:
        - mp = pulsar mass [solar mass]
        - sini = sine of inclination [  ]
        - mf = Keplerian mass function [solar mass]

    Output:
        - mc = companion mass [solar mass]
    """

    mc = 1. # initial guess of 1 solar mass.

    # use one-dimensional Newton-Raphson method.
    for ii in range(n_iter):
        mcb = mc
        f = (mc * sini)**(3./2.) - mf * (mp + mc)
        fp = 3. / 2. * np.sqrt(mc * sini**3) - mf
        mc = mc - f / fp

        if (np.fabs(mc - mcb) < 1e-10):
            break

    return mc

def doppler_shift_orbit(dates, x, pb, ecc, om, t0):
    """
    Computes fraction of period induced by the Dopper effect. 
    (See Equations 8.24, 8.25 of Lorimer & Kramer, 2005.)

    Inputs:
        - x = projected semimajor axis [lt-s]
        - pb = orbital period [days]
        - ecc = eccentricity [  ]
        - om  = argument of periastron [deg]
        - t0 = epoch of periastron passage [MJD] 
        - dates  = epochs to evaluate delay [MJD]

    Output:
        - Doppler factor, v_orb/c [  ]
    """

    ma = mean_anomaly(pb, dates, t0)
    ea = ecc_anomaly(ma, ecc)
    ta = true_anomaly(ea, ecc)
    om = peri_argument(om, pb, ta)
    so, co = np.sin(om * d2r), np.cos(om * d2r)

    amp = 2 * np.pi * x / pb / np.sqrt(1 - ecc**2)
    return amp * (np.cos((om + ta) * d2r) + ecc * co)

def fit_orbit(pars, mjds, periods, periods_errs):
    """
    Applies the Doppler-shifted spin period model to a a set of period measurements and 
    fits for the orbital elements, assuming an eccentric orbit.

    Inputs:
        - pars = a list of the parameters to be fit:
            * pars[0] = spin period [ms]
            * pars[1] = projected semimajor axis [lt-s]
            * pars[2] = orbital period [days]
            * pars[3] = eccentricity [  ]
            * pars[4] = argument of periastron [deg]
            * pars[5] = epoch of periastron passage [MJD]
        - mjds = array of measurement epochs [MJD]
        - periods = measurements of pulsar-spin period [ms]
        - periods_errs = measurement uncertainties [ms]

    Outputs:
        - list of best-fit parameters, same order and units as input parameter list.
    """

    def chisq(parameters):
        p0, p1, p2, p3, p4, p5 = parameters
        model = p0 * (1 + doppler_shift_orbit(mjds, p1, p2, p3, p4, p5))
        return np.sum(((periods - model) / periods_errs)**2)

    result = minimize(chisq, pars, method='Nelder-Mead')
    print "Fit success: {0}".format(result.success)

    pars = result.x
    plt.subplot(211)
    plt.errorbar(mjds, periods, yerr=periods_errs, fmt='ro')
    initial_model = pars[0] * (1 + doppler_shift_orbit(mjds, pars[1], pars[2], pars[3], pars[4], pars[5]))
    plt.plot(mjds, initial_model, 'b-')
    plt.ylabel('Pulse Period (ms)')
    plt.grid()
    plt.subplot(212)
    plt.errorbar(mjds, periods - initial_model, yerr=periods_errs, fmt='ro')
    plt.xlabel('Time (MJD)')
    plt.ylabel('Residual (ms)')
    plt.grid()
    plt.show()

    print "Final parameter estimates:"
    print "    * PS = {0:.6f} ms".format(pars[0])
    print "    * A1 = {0:.6f} lt-s".format(pars[1])
    print "    * PB = {0:.6f} days".format(pars[2])
    print "    * E  = {0:.6f}".format(pars[3])
    print "    * OM = {0:.6f} deg".format(pars[4])
    print "    * T0 = {0:.6f} MJD".format(pars[5])

    return pars

def simOrbit(spin_period, x, pb, ecc, om, t0, n_points=500, n_orbits=10, start_mjd=55720., 
             noise_frac=0.1, uncertainty=1., save=False, outfile="sim_orbit.dat"):

    # simulate time span of data.
    n_days = int(pb * n_orbits)
    mjds = np.linspace(start_mjd, start_mjd + n_days, num=n_points)

    # get time delay due to orbit.

    data = spin_period * (1 + doppler_shift_orbit(mjds, x, pb, ecc, om, t0))
    data += np.random.normal(0., noise_frac * x, n_points)
    uncertainties = np.zeros(n_points) + uncertainty

    if save:

        fout = open(outfile, 'w')
        
        for idx in range(n_points):
            fout.write("{0:.12f}  {1:.12f}  {2:.12f}\n".format(mjds[idx], data[idx], uncertainties[idx]))

        fout.close()

    else:

        return mjds, data, uncertainties
