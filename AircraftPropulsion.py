
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
    in the comments, currently using corrected version
    
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

    """
    gminus = gamma-1
    gminusg = gminus/gamma

    top    = 1 + gminus/2*mach0**2
    bottom = 1 + gminus/2*mach1**2

    first = mach1/mach0*np.sqrt(bottom/top)-1 # Books says -1 should go inside the square root
    first1 = 2*area0area1*first

    second = (top/bottom)**(1/gminusg) - 1
    second1 = 2/(gamma*mach0**2)*second
    
    areamax_area1 = 1 + (first1+second1)/(-pressure_coeff)

    return areamax_area1

def pressure_coeff_crit(mach0, gamma=1.4):
    """Returns the pressure coefficient that would occur when the incoming flow is accelerated sonic conditions"""
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
    gplus = gamma+1
    A0A1a = M1/M0*((1+gminus/2*M0**2)/(1+gminus/2*M1**2))**(gplus/2/gminus)

    return (A0A1a - A0A1, A0A1a)

def main():
    
    arearatios = np.linspace(0, 1, 100)
    arearatio_max = np.empty_like(arearatios)  
    mach0 = 0.9

    cp_crit = pressure_coeff_crit(mach0)
    for i, ar in enumerate(arearatios):
        mach1 = rootfind.newton2(MachfromAreaRatio, .5, M0=mach0, A0A1=ar, gamma=1.4)
        arearatio_max[i] = nacelle_area_ratio(mach0, mach1, ar, cp_crit)

    plt.plot(arearatios, arearatio_max)
    plt.show()

if __name__ == '__main__':
    main()