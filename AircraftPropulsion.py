
import numpy as np
import matplotlib.pyplot as plt
import RootFinding as rootfind
import Component_Design as compdes

def nacelle_area_ratio(mach0, mach1, area0area1, pressure_coeff, gamma:float=1.4):
    """
    Notes
    -----
    This function calculates the area ratio of external cowl area that produces sonic flow on the nacelle to highlight area, if the 
    critical pressure coefficient is given.
    The flow that does not enter the inlet accelerates as it flows around the outside of the cowl and if the cowl is too large the flow
    will become supersonic potentially inducing shock and increasing drag, similar to divergence drag on an airfoil.
    Station numbering 0: freestream, 1: highlight
    p375
    
    Returns
    -------
    0: ratio of maximum-external-area-to-remain-sonic to highlight area

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
    """
    Notes
    -----
    Returns the pressure coefficient that would occur when the incoming flow is accelerated to sonic conditions
    Station numbering 0: freestream, 1: highlight
    p376
    
    Returns
    -------
    pressure coefficient

    Parameters
    ----------
    mach0 : freestream mach number

    Assumptions
    -----------

    """
    gminus = gamma-1
    gplus = gamma+1
    gminusg = gminus/gamma

    top = 1 + gminus/2*mach0**2
    pressure_coeff_crit = 2/(gamma*mach0**2) * ((2*top/gplus)**(1/gminusg) - 1)

    return pressure_coeff_crit

def mach_from_area_ratio(M1:float, M0:float, A0A1:float, gamma:float):
    """
    Notes
    -----
    To be used with newton2 rootfind to determine the mach number to which the flow would be accelerated if it
    experienced an area change from A0 to A1
    Can be reworked to be solved faster and remove dependence on rootfind

    Returns
    -------
    0: difference between input area ratio and calculated area ratio
    1: calculated area ratio

    Parameters
    ----------
    M1 : unknown mach number
    M0 : mach number at known area
    A0A1 : known area ratio
    gamma : ratio of specific heats
    """
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
    p402
    
    Returns
    -------
    ratio of total freestream pressure to fan or compressor face total pressure
    
    Parameters
    ----------
    mach0 : freestream mach number
    """
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
    drag coefficient notes.
    Subscripts i and e denote properties at entry and exit of afterburner.
    Derivation of equation 7.97 for mache^2 is wrong and should be mache^2 = ( A^2ge-2(ge-1)-A\sqrt((A^2+2)ge^2+2) )/(2-A^2)/(ge-1)/ge
    where g is gamma, all other equations are correct for the derivation of Pte/Pti. Figure 7.32 of mache vs machi reflects the 
    incorrect book equation, figure 7.31 cannot be replicated, figure 7.34 cannot be replicated and contradicts figure 7.31.
    The code reflects the correct equation for mach2 presented here
    p506

    Returns
    -------
    ratio of total pressure at the beginning of the afterburner to the total pressure at the exit of the afterburner

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

def prediff_mach(mach2, mach1, Pt31, Ptin, area31, dome_area, gamma):

    first = 1 + (gamma-1)/2*mach1**2
    second = 1 + (gamma-1)/2*mach2**2
    mach22 = Pt31/Ptin*area31/dome_area*mach1*(second/first)**((gamma+1)/(gamma-1)/2)

    return (mach22-mach2,)

def combustor_pressure_loss(Tt31:float, Pt31:float, Tt4:float, dome_vel:float, area31, dome_area, mach31, gamma:float, R_gas:float):
    """
    Notes
    -----
    This function attempts to characterize flow properties in the combustor. There is an empirical relationship for total presure loss due to the prediffuser, heat addition in the combustor, and one for mach number of air to be combusted. Prediffuser pressure loss is from the book, heat comubstion pressure is from class, and mach number is from the book. The class relationship for to-be-combusted mach number for air is in theis function but not used due to uncertainty, but is used in the overall design of the combustor. The book air mach number is higher than the class one and produces a higher change in total pressure. It is recommended to use the book mach number.
    Station numbering "in" entrance of dome
    p517

    Returns
    -------
    0: the change in total pressure from the pre combusted air to combusted gas mixture (Pt_combusted - Pt_not comubsted)
    1: total pressure of air entering combustor
    2: mach number of air to be combusted

    Parameters
    ----------
    Tt31 : 
    Pt31 : 
    Tt4 : 
    dome_vel : velocity in the dome
    gamma : ratio of specific heats of the air
    R_gas : gas constant of the air

    Assumptions 
    -----------
    0: prediffuser design in a dump diffuser
    1: an area ratio 1<A2/A1<5 
    """

    first = (1-area31/dome_area)
    Ptin = Pt31*np.exp(-gamma*mach31**2/2*(first**2+first**6))
    
    # machin = dome_vel/np.sqrt(gamma*R_gas*Tt31)
    deltaPt = Ptin*0.53*machin**2*(0.95+.05*(Tt4/Tt31))
    print(deltaPt, machin)

    machin_book = rootfind.newton2(prediff_mach, .5, mach1=mach31, Pt31=Pt31, Ptin=Ptin, area31=area31, dome_area=dome_area, gamma=gamma) # derived from book p517
    deltaPt = Ptin*0.53*machin_book**2*(0.95+.05*(Tt4/Tt31))
    print(deltaPt, machin_book)

    return deltaPt, Ptin, machin_book, machin

def main():
    combustor_values = compdes.combustor_design_test_values()
    compdes.compressor_design_test_values()
    ref_height = combustor_values[-4]
    dome_vel = combustor_values[-3]
    
    mach31 = .3
    Tt31 = 897.1*5/9 # R
    Pt31 = 52.64*6895 # psia
    pitch_diam = 28.02*.0254# in
    # Combustor parameters
    Tt4 = 2560*5/9 # R
    area31 = 0.08868156921060048
    dome_area = ref_height*np.pi*pitch_diam

    dPt, Ptin, machin2, machin = combustor_pressure_loss(Tt31, 1, Tt4, dome_vel, area31, dome_area, .5, 1.4, 287.05)

if __name__ == '__main__':
    main()