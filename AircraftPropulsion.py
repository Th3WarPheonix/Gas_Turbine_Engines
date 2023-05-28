
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

    p375
    
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

    first = mach1/mach0*np.sqrt(firsta/firstb)-1
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

def mach_from_area_ratio(M1, M0, A0A1, gamma):
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

def supersonic_pressure_recovery(mach0, std='MIL'):
    """
    Notes
    -----
    Pressure recovery of supersonic inlets according to two empirical models 1: military standard Mil-E_5008B and
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

def afterburner_pressure_loss(machi, drag_coeff, gammai, gammae, q=None, cpi=None, Ti=None, mode='dry'):
    """
    Notes
    -----
    Total pressure loss in a constant area afterburner in either dry or wet afterburner
    drag coefficient notes
    subscripts i and e denote properties at entry and exit of afterburner
    Derivation of equation 7.97 for mache^2 is wrong and should be mache^2 = ( A^2ge-2(ge-1)-A\sqrt((A^2+2)ge^2+2) )/(2-A^2)/(ge-1)/ge
    where g is gamma, all other equations are correct for the derivation of Pte/Pti. Figure 7.32 of mache vs machi reflects the 
    incorrect book equation, figure 7.31 cannot be replicated, figure 7.34 cannot be replicated and contradicts figure 7.31.
    The code reflects the correct equation for mach2 presented here
    p506

    Parameters
    ----------
    machi : mach number at beginning of afterburner
    drag_coeff : drag coefficient of the flame holders
    gammai : ratio of specific heats at beginning of afterburner
    gammae : ratio of specific heats at exit of afterburner
    dry : no fuel combusted in afterburner
    wet : fuel combusted in afterburner
    q : amount of heat released when fuel is combusted
    cpi : specific heat at constant pressure of exhuast entering afterburner
    Ti : temperature of exhaust entering afterburner

    Assumptions
    -----------
    0: neglects fuel-air-ratio
    """
    giminus = gammai-1
    geminus = gammae-1
    geminusg = geminus/gammae
    giminusg = giminus/gammai

    fourth = (1+gammai*machi**2*(1-drag_coeff/2))/gammai/machi
    fifth = giminus/(1+giminus/2*machi**2)
    if mode=='wet': eight = 1+q/cpi/Ti
    else: eight = 1
    third = fourth*np.sqrt(fifth/eight)

    second = (2-third**2)*geminus*gammae
    mache_sq = ( third**2*gammae-2*geminus-third*np.sqrt((third**2-2)*gammae**2+2) )/second

    sixth = 1+gammai*machi**2*(1-drag_coeff/2)
    seventh = 1+gammae*mache_sq
    pressure_ratio = sixth/seventh*(1+geminus/2*mache_sq)**(1/geminusg)/(1+giminus/2*machi**2)**(1/giminusg)

    return pressure_ratio

def combustor_pressure_loss():
    """Notes
    empirical relationship for total presure loss due to heat addition in a combustor"""

    deltaPt = Ptin*0.53*machin**2*(0.95+.05*(Tt4/Tt31))

    return deltaPt

def main():
    machis = np.linspace(0.1, .5, 1000)
    prs = np.empty_like(machis)
    CD = 1.25
    # print(afterburner_pressure_loss(0.25, CD, 1.33, 1.33, 1, 1, 1))
    for i, machi in enumerate(machis):
        prs[i] = afterburner_pressure_loss(machi, CD, 1.33, 1.33, 1, 1, 1)

    plt.plot(machis, prs)
    plt.show()

if __name__ == '__main__':
    main()