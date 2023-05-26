
import numpy as np
import matplotlib.pyplot as plt
import RootFinding as rootfind

def nacelle_area_ratio(mach0, mach1, area0area1, pressure_coeff, gamma:float=1.4):
    """
    Notes
    -----
    This function calculates the area ratio of external cowl area that produces sonic flow on the nacelle to highlight area, if the 
    critical pressure coefficient is given.
    The flow that does not enter the inlet accelerates as it flows around the outside of the cowl and if the cowl is too large the flow
    will become supersonic potentially inducing shock and increasing drag, similar to divergence drag on an airfoil.

    Derivation on pages 375 and 376 of Aircraft Propulsion 3e is wrong for this equation in the spot noted 
    in the comments, currently using incorrect version to verify rest of calculations
    
    Returns
    -------
    ratio of maximum external area to remain sonic to highlight area

    Parameters
    ----------
    mach0 : freestream mach number
    mach1 : mach number at highlight
    area0area1 : ratio of streamtube area to highlight area
    pressure_coeff : pressure coefficient that occurs when flow is sonic

    Assmuptions
    -----------
    0: compressible flow
    """
    gminus  = gamma-1
    gminusg = gminus/gamma

    firsta = 1 + gminus/2*mach0**2
    firstb = 1 + gminus/2*mach1**2

    first = mach1/mach0*np.sqrt(firsta/firstb)-1 # Books says firsta should be in the denominator
    first1 = 2*area0area1*first

    second = (firsta/firstb)**(1/gminusg) - 1
    second1 = 2/(gamma*mach0**2)*second
    
    areamax_area1 = 1 + (first1+second1)/(-pressure_coeff)

    return areamax_area1

def pressure_coeff_crit(mach0, gamma=1.4):
    """Returns the pressure coefficient that would occur when the incoming flow is accelerated sonic conditions
    0 freestream
    1 highlight
    p376"""
    gminus = gamma-1
    gplus = gamma+1
    gminusg = gminus/gamma

    top = 1 + gminus/2*mach0**2
    pressure_coeff_crit = 2/(gamma*mach0**2) * ((2*top/gplus)**(1/gminusg) - 1)

    return pressure_coeff_crit

def MachfromAreaRatio(M1, M0, A0A1, gamma):
    """To be used with newton2 rootfind to determine the mach number to which the flow would be accelerated if it
    expereinced an area change from A0 to A1"""
    gminus = gamma-1
    gplus  = gamma+1
    A0A1a = M1/M0*((1+gminus/2*M0**2)/(1+gminus/2*M1**2))**(gplus/2/gminus)

    return (A0A1a - A0A1, A0A1a)

def additive_drag(mach0, mach1, gamma=1.4):
    """Additive drag due to air intake
    p378"""
    gminus  = gamma-1
    gplus   = gamma+1
    gminusg = gminus/gamma

    first0 = 1 + gminus/2*mach0**2
    first1 = 1 + gminus/2*mach1**2
    first = (first0/first1)**(gplus/2/gminus)

    second = mach1*np.sqrt(first0/first1)-mach0

    third = (first0/first1)**(1/gminusg)

    drag = gamma*mach1*first*second+third-1

    return drag

def pressure_recovery_supersonic_inlet(mach0, std='MIL'):
    """
    Notes
    -----
    Pressure recovery of supersonic inlets according to 1: military standard Mil-E_5008B and
    2: AIA (Aircraft industries Association)
    Both standards are considered conservative (worst case pressure recovery) by todays standards

    p402"""
    match std:
        case 'MIL':
            if 1 <= mach0 < 5: pressure_ratio = 1-.075*(mach0-1)**1.35
            else: pressure_ratio = 800/(mach0**4+935)
        case 'AIA':
            pressure_ratio = 1-.1*(mach0-1)**1.5

    return pressure_ratio

def main():
    mach0s = np.linspace(1, 5, 100)
    prmil = np.empty_like(mach0s)
    praia = np.empty_like(mach0s)
    for i, m0 in enumerate(mach0s):
        prmil[i] = pressure_recovery_supersonic_inlet(m0)
        praia[i] = pressure_recovery_supersonic_inlet(m0, 'AIA')

    plt.plot(mach0s, prmil, label='mil')
    plt.plot(mach0s, praia, label='aia')
    plt.show()

if __name__ == '__main__':
    main()