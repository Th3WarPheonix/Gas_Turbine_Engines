
"""Module containing basic isentropic flow relation equations e.g. stagnation temrperature from Mach number and pressure and vice versa, for all properties"""

def total_pressure(M, p, gamma=1.4):
    """Returns stagnation pressure from static pressure"""
    P0P1 = (1 + (gamma-1)/2.0*M**2.0)**(gamma/(gamma-1))
    return P0P1*p

def total_temperature(M, T, gamma=1.4):
    """Returns stagnation temperature from static temperature"""
    T0T1 = (1 + (gamma-1)/2.0*M**2.0)
    return T0T1*T

def total_density(M, r, gamma=1.4):
    """Returns stagnation density from static density"""
    rho0rho1 = (1 + (gamma-1)/2.0*M**2.0)**(1.0/(gamma-1))
    return rho0rho1*r

def static_pressure(M, p0, gamma=1.4):
    """Returns static pressure from stagnation pressure"""
    P0P1 = (1 + (gamma-1)/2.0*M**2.0)**(gamma/(gamma-1))
    return P0P1**-1*p0

def static_temperature(M, T0, gamma=1.4):
    """Returns static temperature from stagnation temperature"""
    T0T1 = (1 + (gamma-1)/2.0*M**2.0)
    return T0T1**-1*T0

def static_density(M, r0, gamma=1.4):
    """Returns static density from stagnation density"""
    rho0rho1 = (1 + (gamma-1)/2.0*M**2.0)**(1.0/(gamma-1))
    return rho0rho1**-1*r0

def pressure_ratio(M, gamma=1.4):
    """Returns pressure ratio, P1/P0"""
    P0P1 = (1 + (gamma-1)/2.0*M**2.0)**(gamma/(gamma-1))
    return P0P1**-1

def temperature_ratio(M, gamma=1.4):
    """Returns temperature ratio, T1/T0"""
    T0T1 = (1 + (gamma-1)/2.0*M**2.0)
    return T0T1**-1

def density_ratio(M, gamma=1.4):
    """Returns stagnation density ratio, rho1/rho0"""
    rho0rho1 = (1 + (gamma-1)/2.0*M**2.0)**(1.0/(gamma-1))
    return rho0rho1**-1

def pressure_ratio2mach(ratio, gamma=1.4):
    """Ratio is given as static pressure/total pressure"""
    return ((ratio**(-(gamma-1)/gamma)-1)*2/(gamma-1))**0.5

def temperature_ratio2mach(ratio, gamma=1.4):
    """Ratio is given as static temperature/total temperature"""
    return ((1/ratio-1)*2/(gamma-1))**0.5

def density_ratio2mach(ratio, gamma=1.4):
    """Ratio is given as static density/total density"""
    return ((ratio**(-(gamma-1)/1)-1)*2/(gamma-1))**0.5