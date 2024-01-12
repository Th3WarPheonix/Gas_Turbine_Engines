
import numpy as np
import USAtmos1976
import Unit_Conversions as units
import matplotlib.pyplot as plt

Pstd = 101325 # Pa
Tstd = 288.15 # K
densitystd = 1.225 # kg/m^3
spdsnd_std = 340.3 # m/s

"""
Variables 
---------
variable name : actual name : quantity/description : in book

wngld : wing loading : takeoff_weight/wing_planform_area : Wto/S
thstld : thrust laoding : thrust/takeoff_weight : T/Wto
engtype : engine type : tj, lbtf, tp, etc : 
engmode : engine mode : wet or dry : 
thrtlrto : throttle ratio :  : TR
thstlps : installed thrust lapse : thrust/sea_level_thrust: thstlps
instwf : instantaneous weight fraction : weight/takeoff_weight : instwf

alt : altiude : : 
gamma : ratio of specific heat : cp/cv : gamma
Rgas : gas specific gas constant : : R

CLmax : maximum coefficient of lift : lift/density/S/V^2/2 : C_Lmax
CD : coefficient of drag : drag/density/S/V^2/2 : C_D
CD0 : parasite drag : drag coeff at zero lift, CDmin + K"Cmin^2: C_D0
CDR : coefficient of additional drags : : C_DR
CDstar : drag coefficient at max L/D : : C*_D
Kp1 : invscid/induced drag coefficient : 1/pi/aspct_rto/efficiency : K'
Kp2 : viscid drag coefficient : K" : K"
K1 : drag polar coefficient 1 : K1 = K' + K" : K1
K2 : drag polar coefficient 2 : K2 = -2K"CLmin : K2

drgcf : drag coefficient :  : epsilon
frccf : friction coefficient :  : mu
velrt : velocity ratio : velocity/stall_velocity : k

thetat : nondimensional total temperature : Tt/Tstd : theta_0
thetas : nondimensional static temperature : Ts/Tstd : theta
deltat : nondimensional total pressure : Pt/Pstd : delta_0
deltas : nondimensional static pressure : Ps/Pstd : delta
sigmat : nondimensional total density : densityt/densitystd : sigma_0
sigmat : nondimensional static density : densitys/densitystd : sigma
dynpress : dynamic pressure : density*velocity^2/2 : q
thetabrk : theta break : control system maximum Tt4 and compressor 
    pressure ratio : theta_0_break

swpang : wing sweep angle :  : Lambda
ldfactor : load factor :  : n
thstang : angle of thrust vector to chord line :  : psi 
emptwf : empty aircraft weight fraction : We/Wto : Gamma
drgthst : total drag-to-thrust ratio : : u
tsfc : thrust specific fuel consumption : fuel_rate/thrust : tsfc
mfp : mass flow parameter : fcn(gas, mach) : mfp 

Abbreviations 
-------------
BCA : best cruise altitude
BCM : best cruise altitude
SLS : sea level static

Subscripts
----------
dry : afterburner off
wet : afterburner on
anlys : analysis
horz : horizontal
accel : acceleration
ln : landing
to : takeoff
td : touchdown
obs : obstacle
drgplr : drag polar

"""

def get_mfp(mach, gamma, Rgas):
    """
    Returns
    -------
    mfp : mass flow parameter m*sqrt(Tt)/Pt/area
    """
    exponent = (gamma+1)/2/(1-gamma)
    mfp = mach*np.sqrt(gamma/Rgas)*((1+(gamma-1)/2)*mach**2)**exponent
    return mfp

def get_atmos(altitude, day='std'):
    """
    Notes
    -----
    Altitude in feet. Can be a numpy array
    
    Returns
    -------
    0: thetas Ts/Tstd
    1: deltas Ps/Pstd
    2: temperature
    3: pressure
    4: density
    """
    match day:

        case 'std':
            temp = USAtmos1976.stdday(altitude*.0254*12).temperature[0]
            press = USAtmos1976.stdday(altitude*.0254*12).pressure[0]
            density = USAtmos1976.stdday(altitude*.0254*12).pressure[0]
        case 'hot':
            temp = USAtmos1976.hotday(altitude*.0254*12).temperature[0]
            press = USAtmos1976.hotday(altitude*.0254*12).pressure[0]
            density = USAtmos1976.hotday(altitude*.0254*12).pressure[0]
        case 'trop':
            temp = USAtmos1976.tropday(altitude*.0254*12).temperature[0]
            press = USAtmos1976.tropday(altitude*.0254*12).pressure[0]
            density = USAtmos1976.tropday(altitude*.0254*12).pressure[0]
        case 'cold':
            temp = USAtmos1976.coldday(altitude*.0254*12).temperature[0]
            press = USAtmos1976.coldday(altitude*.0254*12).pressure[0]
            density = USAtmos1976.coldday(altitude*.0254*12).pressure[0]
        case _:
            print('not a known type of day')
            return None
    
    thetas = temp/288.15
    deltas = press/101325
    return (thetas, deltas, temp, press, density)

def get_dynpress(alt, mach, gamma=1.4):
    """
    Notes
    -----
    Altitude in feet. Can be a nummpy array
    
    pg 11
    
    Returns
    -------
    dynamic pressure in Pascals

    Parameters
    ----------
    alt : altitude in feet
    """

    deltas = get_atmos(alt, 'std')[1]
    return deltas*101325*mach**2*gamma/2


def get_drgplr(mach0:int, CLmax, velrt=1.25, K2=0, CL=0, suppress=False):
    """
    Notes
    -----
    Calculate drag polar
    K1, K2, CD0 are assumed known
    high performance aircraft CLmin~=0 K2~=0
    K2 = 0 for uncambered
    .001 <= K" <= .03
    .1 <= CLmin <= .3
    7 <= AR <= 10
    .75 <= efficiency <= .85

    Lift-drag polar breakdown
    CD = CDmin + K'CL^2 + K"(CL - CLmin)^2
    CD = K1CL^2 + K2CL + CDO
    
    pg 25,37

    Returns
    ------
    0: CD0
    1: K1
    2: CD

    Parameters
    ----------
    mach0 : freestream mach number
    CLmax : maximum coefficient of lift
    velrt : velocity ratio V = velrt*Vstall
    """

    if 0 <= mach0 <= 1:
        K1 = 0.18
    elif 1 < mach0 <= 2:
        K1 = (.18 - .36)/(1-2)*(mach0-1)+.18
    else:
        K1 = 0
        if not suppress: print("Flight mach number not in range for turbine engines1")

    if 0 <= mach0 <= 0.8:
        CD0 = .014
    elif  0.8 < mach0 <= 0.9:
        CD0 = (.014 - .016)/(0.8-.9)*(mach0-0.8)+.014
    elif 0.9 < mach0 <= 1.1:
        CD0 = (.016 - .026)/(0.9-1.1)*(mach0-0.9)+.016
    elif 1.1 < mach0 <= 1.2:
        CD0 = (.026 - .028)/(1.1-1.2)*(mach0-1.1)+.026
    elif 1.2 < mach0 <= 2:
        CD0 = .028
    else:
        CD0 = 0
        if not suppress: print("Flight mach number not in range for turbine engines2")
        
    if not CL: CL = CLmax/velrt**2
    CD = K1*CL**2+K2*CL+CD0

    return (CD0, K1, CD)

def get_thstlps(thrtlrto, mach0, engine_type, engmode, day, alt, gamma=1.4):
    """
    Notes
    -----
    thstlps is the installed full throttle thrust lapse, which depends on
    altitude, speed, and whether or not an afterburner is operating. i.e
    thstlps = thrust/thrust_sea_level

    pg 38
    
    Parameters
    ----------
    engine_type : 
        hbtf : High bypass turbofan
        lbtf : Low bypass turbofan
        tj : Turbojet
        tp : Turboprop
    engmode : 
        mil : military power
        max : max power
    alt : altitude in feet
    """

    thetas = get_atmos(alt, day)[0]
    deltas = get_atmos(alt, day)[1]

    # normalized total temperature and total pressure
    norm_Tt = thetas*(1+(gamma-1)/2*mach0**2)
    norm_Pt = deltas*(1+(gamma-1)/2*mach0**2)**(gamma/(gamma-1))

    if norm_Tt <= thrtlrto:
        zeta = 0
    else:
        zeta = 1

    match engine_type:
        case 'hbtf': # high bypass turbofan
            thstlps = norm_Pt*(1-0.49*np.sqrt(mach0) - 3*(norm_Tt-thrtlrto)/
                                    (1.5+mach0)*zeta)
        case 'lbtf': # low bypass turbofan
            if engmode == 'mil': # military power
                thstlps = 0.6*norm_Pt*(1-3.8*(norm_Tt-thrtlrto)/norm_Tt*zeta)
            elif engmode == 'max': # maximum power
                thstlps = norm_Pt*(1-3.5*(norm_Tt-thrtlrto)/norm_Tt*zeta)
        case 'tj': # turbojet
            if engmode == 'mil': # military power
                thstlps = 0.8*norm_Pt*(1-.16*np.sqrt(mach0) - 
                                            25/norm_Tt*(norm_Tt-thrtlrto)
                                            /(9+mach0)*zeta)
            elif engmode == 'max': # maximum power
                thstlps = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(mach0) -
                                         1.5*(norm_Tt-thrtlrto)/norm_Tt*zeta)
        case 'tp': # turboprop
            thstlps = norm_Pt # less than mach = 0.1
            thstlps = norm_Pt*(1-0.96*(mach0-1)**0.25 - 3/8.13*(
                norm_Tt-thrtlrto)/(mach0-0.1)*zeta)
    
    return thstlps


def get_takeoff_vel(instwf, CLmax, densitys, wngld, kto=1.2):
    """
    Notes
    -----
    Calculates the takeoff velocity
    
    p 29
    """
    vel_to = kto*np.sqrt(2*instwf/densitys/CLmax*wngld)
    print(284, kto, instwf, densitys, CLmax, wngld)
    return vel_to


def thstld_takeoff(wngld, densitys, CLmax, mach, thstlps, instwf, frccf_ln, 
                   velrt_to, takeoff_distance, rot_time=3, CL=0, CDR=0, 
                   g0=9.81):
    """
    Notes
    -----
    Thrust loading during takeoff acceleration up to rotation
    instantaneous weight is
    given by
    W = instwf*WTO
    instwf, instantaneous weight fraction, depends on how much fuel has
    been consumed and payload delivered; more AED2e p24

    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    standard rotation time is 3 secs
    epsilon is the sum of the drags on the aircraft
    kto is the ratio of takeoff veloicty to stall velocity

    pg 34
    """
    CD = get_drgplr(mach, CLmax, velrt_to)[2]

    frccf_to = CD + CDR + frccf_ln*CL
    help3 = instwf/densitys/g0/frccf_to
    help2 = rot_time*velrt_to*np.sqrt(2*instwf/densitys/CLmax)
    help1 = -(takeoff_distance-help2*np.sqrt(wngld))/help3/wngld
    thstld = (frccf_to*velrt_to**2/(1-np.exp(help1))/CLmax+
                      frccf_ln)*instwf/thstlps

    return thstld

def thstld_cruise(wngld, alt, mach, CLmax, thstlps, instwf, 
                       CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during cruise
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 26
    """
    dyn_press = get_dynpress(alt, mach)
    CD0, K1, CD = get_drgplr(mach, CLmax)
    thstld = instwf/thstlps*(K1*instwf/dyn_press*wngld + K2 + (CD0+CDR)/instwf*
                         dyn_press/wngld) 

    return thstld

def thstld_turn(wngld, load_factor, alt, mach, CLmax, thstlps, instwf, CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 26
    """
    dyn_press = get_dynpress(alt, mach)
    CD0, K1, CD = get_drgplr(mach, CLmax)

    help1  = K1*load_factor**2*instwf/dyn_press
    help2 = (CD0+CDR)/instwf*dyn_press
    thstld = instwf/thstlps*(help1*wngld + K2*load_factor + 
                                 help2/wngld) 

    return thstld

def thstld_horzaccel(wngld, deltamach, deltat, air_temp, alt, mach, CLmax, 
                     thstlps, instwf, CDR=0, K2=0, g0=9.81, gamma=1.4, Rgas=287.05):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 27
    """
    dyn_press = get_dynpress(alt, mach)
    CD0, K1, CD = get_drgplr(mach, CLmax)

    help1  = K1*instwf/dyn_press
    help2 = (CD0+CDR)/instwf*dyn_press
    help3 = np.sqrt(gamma*Rgas*air_temp)*(deltamach)/g0/deltat
    thstld = instwf/thstlps*(help1*wngld + K2 + 
                                 help2/wngld + help3) 

    return thstld

def thstld_landing(wngld, densitys, CLmax, CD, thstlps, instwf, frccf_ln, velrt, distance, rotation_time=3, 
                        CL=0, CDR=0, g0=9.81):
    """
    Notes
    -----
    Thrust loading during landing free roll and braking
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    standard rotation time is 3 secs
    epsilon is the sum of the drags on the aircraft
    kto is the ratio of takeoff veloicty to stall velocity

    pg 31
    """
    epsilon = CD + CDR + frccf_ln*CL
    help3 = instwf/densitys/g0/epsilon
    help2 = rotation_time*velrt*np.sqrt(2*instwf/densitys/CLmax)
    
    if thstlps != 0:
        help1 = (distance-help2*np.sqrt(wngld))/help3/wngld
        thstld = (epsilon*velrt**2/(np.exp(help1)-1)/CLmax-
                          frccf_ln)*instwf/thstlps
    else:
        fourth = help3*np.log(1+epsilon/frccf_ln/CLmax*velrt**2)
        thstld = ((-help2+np.sqrt(help2**2+4*fourth*distance))/
                          2/fourth)**2

    return thstld

def excess_power(wngld, thstld, load_factor, velocity, alt, mach, CLmax, thstlps,
                 instwf, K2=0, CDR=0):
    """
    Notes
    -----
    Calculates excess power for the aircraft at the defined thrust
    loading, wing loading, throttle ratio for any instwf, load factor
    across the flight envelope (altitude vs altitude)

    pg 48

    Returns
    -------
    excess power

    Parameters
    ----------
    alt : altitude in feet no regard for other units
    all other parameters units in SI
    """
    dyn_press = get_dynpress(alt, mach)
    CD0, K1, CD = get_drgplr(mach, CLmax, suppress=True)

    help1  = thstlps/instwf*thstld
    help2 = -K1*load_factor**2*instwf/dyn_press*wngld
    help3  = -K2*load_factor
    fourth = -(CD0+CDR)/instwf/wngld*dyn_press
    exs_pwr = velocity*(help1+help2+help3+fourth)

    return exs_pwr 

def constraint_anlys(engtype):
    """
    Notes
    -----
    Need to aid remining constraint equations form chapter 2 but
    currently not necessary for remainder of book walkthrough
    """
    
    wngld1a = np.linspace(20, 120, 30)
    wngld1 = units.convert_pressure(wngld1a/144, 'Pa')
    
    densitys = 1.05498044 # kg/m^3
    TR = 1.07
    instwf = [1, .78, .78, .78, .78, .56, .78]
    takeoff_distance = 1500*12*.0254 # m
    mach = [.1, 1.5, 1.6, .9, 1.2, 0, 1.8]
    alt = [2000, 30000, 30000, 30000, 30000, 0, 40000]

    load_turn1 = load_turn2 = 5 # g

    deltamach = .8
    deltat = 50

    air_temp = get_atmos(30000, 'std')[0]*288.15
    kto = 1.2
    ktd = 1.15
    muto = .05
    mutd = .18
    CLmax = 2

    landing_drag = .8123
    thrust_reverse = 0

    thstlps1 = get_thstlps(TR, mach[0], engtype, 'max', alt[0])
    thstlps2 = get_thstlps(TR, mach[1], engtype, 'mil', alt[1])
    thstlps3 = get_thstlps(TR, mach[2], engtype, 'max', alt[2])
    thstlps4 = get_thstlps(TR, mach[3], engtype, 'max', alt[3])
    thstlps5 = get_thstlps(TR, mach[4], engtype, 'max', alt[4])
    thstlps7 = get_thstlps(TR, mach[6], engtype, 'max', alt[6])
    # All thrust loading equations have been validated
    thstld1 = thstld_takeoff(wngld1, densitys, CLmax, mach[0], thstlps1, 
                                 instwf[0], muto, kto, takeoff_distance) 
    thstld2 = thstld_cruise(wngld1, alt[1], mach[1], CLmax, thstlps2,
                                instwf[1]) 
    thstld3 = thstld_turn(wngld1, load_turn1, alt[2], mach[2], 
                              CLmax, thstlps3, instwf[2]) 
    thstld4 = thstld_turn(wngld1, load_turn2, alt[3], mach[3], CLmax, 
                              thstlps4, instwf[3]) 
    thstld5 = thstld_horzaccel(wngld1, deltamach, deltat, air_temp,
                                    alt[4], mach[4], CLmax, thstlps5, instwf[4]) 
    thstld6 = thstld_landing(wngld1, densitys, CLmax, landing_drag, 
                                 thrust_reverse, instwf[5], mutd, ktd, 
                                 takeoff_distance)
    thstld7 = thstld_cruise(wngld1, alt[6], mach[6], CLmax, thstlps7, 
                                instwf[6]) 
    
    thstld6 = units.convert_pressure(thstld6, "psi")*144

    plt.plot(wngld1a, thstld1, label="Takeoff", linestyle='--')
    plt.plot(wngld1a, thstld2, label="Cruise", linestyle='--')
    plt.plot(wngld1a, thstld3, label="Turn 1", linestyle='--')
    plt.plot(wngld1a, thstld4, label="Turn 2", linestyle='--')
    plt.plot(wngld1a, thstld5, label="Accel", linestyle='--')
    plt.vlines(thstld6, .4, 1.6, label="Landing", linestyle='--')
    plt.plot(wngld1a, thstld7, label="Max Mach", linestyle='--')

    plt.legend()
    plt.ylim([.4, 1.6])
    plt.xlim([20, 120])
    plt.show()
    plt.close()

def power_anlys(thstld, wngld, thrtlrto, engtype, CLmax):
    """
    Notes
    -----
    Function to compare this codes results to book results
    
    Mimics plot on pg 48 but slight indent on righgt corner of envelope
    Single point calculation performed after plotting gives 316 ft/s
    book gives 320 ft/s

    More resolution needed but takes to long to render plot, convert to
    C or Fortran
    """

    load_factor = 1
    instwf = .97
    wngld = units.convert_pressure(wngld/144, 'Pa')
    N = 50
    velocity = np.linspace(100, 1900, N)*.0254*12 # ft/sec
    altitude = np.linspace(0, 60000, N) # ft
    power = np.zeros((N, N))
     
    for i, alt in enumerate(altitude):
        for j, vel in enumerate(velocity):
            mach = (vel/USAtmos1976.stdday(alt).speed_of_sound)[0]
            thstlps = get_thstlps(thrtlrto, mach, engtype, 'mil', alt)
            power[i][j] = excess_power(wngld, thstld, load_factor, vel, alt, 
                                       mach, CLmax, thstlps, instwf)
    power = power/12/.0254

    levels = np.arange(0, 600+1, 50)
    plt.contourf(velocity/12/.0254, altitude, power, levels)
    plt.colorbar()
    plt.xlabel('Velocity (ft/s)')
    plt.ylabel('Alt (ft)')
    plt.savefig('Excess Power Envelope.png')

    alt = 36000
    vel = 1400*12*.0254
    mach = vel/(atmos(alt).speed_of_sound)[0]
    thstlps = get_thstlps(thrtlrto, mach, engtype, 'mil', alt)
    pwr = excess_power(wngld, thstld, load_factor, vel, alt, mach, 
                       CLmax, thstlps, instwf)
    print(pwr/.0254/12)

def get_empwf(takeoff_weight, type='fighter'):
    """
    Notes
    -----
    Preliminary results for empty weight fraction, the ratio of empty
    aircrat weight to takeoff weight (We/Wto).

    p 71
    
    Parameters
    ---------
    type : type of aircraft
        fighter cargo passenger twin tprop
    """
    match type:
        case 'fighter':
            empwf = 2.34*takeoff_weight**-0.13
        case 'cargo':
            empwf = 1.26*takeoff_weight**-0.08
        case 'passenger':
            empwf = 1.02*takeoff_weight**-0.06
        case 'twin tprop':
            empwf = 0.96*takeoff_weight**-0.05

    return empwf

def tsfc_initial(mach0, theta, engtype, engmode):
    """
    Notes
    -----
    Provide an initial estimate for the thrust specific fuel consumption
    units of 1/hour

    p 71
    """
    match engtype:
        case 'hbtf': # high bypass turbofan
            tsfc = (0.45+0.54*mach0)*np.sqrt(theta)
        case 'lbtf': # low bypass turbofan
            if engmode == 'mil': # military power
                tsfc = (0.9+0.30*mach0)*np.sqrt(theta)
            elif engmode == 'max': # maximum power
                tsfc = (1.6+0.27*mach0)*np.sqrt(theta)
        case 'tj': # turbojet
            if engmode == 'mil': # military power
                tsfc = (1.1+0.30*mach0)*np.sqrt(theta)
            elif engmode == 'max': # maximum power
                tsfc = (1.5+0.23*mach0)*np.sqrt(theta)
        case 'tp': # turboprop
            tsfc = (0.18+0.8*mach0)*np.sqrt(theta)

    return tsfc

def get_engconst(engtype, engmode):
    """
    Notes
    -----
    Get engine type specific constants for mission analysis

    p 71

    Parameters
    ----------
    engtype:
        hbtf : high bypass turbofan constants returned
        lbtf : high bypass turbofan constants returned
        tj : turbojet constants returned
        tp : turboprop constants returned
    engmode:
        mil : military power constants returned
        max : maximum power constants returned

    Returns
    -------
    0: C1 units 1/hour
    1: C2 units 1/hour

    """
    match engtype:
        case 'hbtf': # high bypass turbofan
            consts = [0.45, 0.54]
        case 'lbtf': # low bypass turbofan
            if engmode == 'mil': # military power
                consts = [0.9, 0.30]
            elif engmode == 'max': # maximum power
                consts = [1.6, 0.27]
        case 'tj': # turbojet
            if engmode == 'mil': # military power
                consts = [1.1, 0.30]
            elif engmode == 'max': # maximum power
                consts = [1.5, 0.23]
        case 'tp': # turboprop
            consts = [0.18, 0.8]
    
    return consts

def get_bcmbca(instwf, CLmax, wngld, machcrit=0.9, CDR=0, Pstd=101325, gamma=1.4):
    """
    Notes
    -----
    Obtain the best cruising altitude for subsonic speeds. The best
    cruise mach number is M=.9 after which divergence drag begins
    """
    CD0, K1, _ = get_drgplr(machcrit, CLmax)
    bca = 2*instwf*wngld/gamma/Pstd/machcrit*2*np.sqrt((CD0+CDR)/K1)

    return bca

def wf_warmup(alt, deltat, thstld, engtype, engmode, thstlps, instwf):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after warm-up

    p 67
    """
    thetas = get_atmos(alt)[0]
    C1, _ = get_engconst(engtype, engmode)

    wght_frac = 1 - C1/3600*np.sqrt(thetas)*thstlps/instwf*thstld*deltat
    
    return wght_frac

def wf_toaccel(alt, drgcf, frccf, mach, mach_to, vel, thstld, wngld, engtype,
                      engmode, thstlps, instwf, g0=9.81):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff acceleration

    p 63
    """
    thetas = get_atmos(alt)[0]
    C1, C2 = get_engconst(engtype, engmode)
    dyn_press = get_dynpress(alt, mach_to)

    ttldrag_thrst = (drgcf*dyn_press/instwf/wngld+frccf)*(instwf/thstlps/
                                                          thstld)
    print(700, instwf, thstlps, thstld)
    wght_frac = np.exp(-(C1+C2*mach)/3600*np.sqrt(thetas)/g0*vel/
                        (1-ttldrag_thrst))

    return wght_frac

def wf_torot(alt, mach_to, thstld, engtype, engmode, thstlps, instwf, 
                    rot_time=3):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff rotation

    p 68
    """
    thetas = get_atmos(alt)[0]
    C1, C2 = get_engconst(engtype, engmode)

    wght_frac = 1-(C1+C2*mach_to)/3600*np.sqrt(thetas)*thstlps/instwf*(
        thstld*rot_time)

    return wght_frac

def wf_takeoff(alt_to, thstld, wngld, thrtlrto, engtype, CLmax):
    """
    Notes
    -----
    Calculates the weight fraction the takeoff sequence: warm-up, acceleration, and rotation
    Returns three weight fraction values: warm-up, acceleration, and rotation
    """
    
    mach_wu = 0
    instwf_wu = 1
    deltat = 60 # s
    engmode_wu = 'mil'
    thstlps_wu = get_thstlps(thrtlrto, mach_wu, engtype, engmode_wu, 'std', alt_to)
    instwf_wu = wf_warmup(alt_to, deltat, thstld, engtype, engmode_wu, thstlps_wu, instwf_wu)

    mach_to = .1819
    mach_toa = .1
    drgcf_toa = .3612
    frccf_toa = .05
    engmode_toa = 'max'
    densitys = get_atmos(alt_to, 'std')[4]
    vel_to = get_takeoff_vel(instwf_wu, CLmax, densitys, wngld)
    thstlps_toa = get_thstlps(thrtlrto, mach_toa, engtype, engmode_toa, 'std', alt_to)
    instwf_toa = wf_toaccel(alt_to, drgcf_toa, frccf_toa, mach_toa, mach_to, vel_to, thstld, wngld, engtype,
                      engmode_toa, thstlps_toa, instwf_wu)
    
    instwf1 = instwf_toa*instwf_wu
    thstlps_torot = get_thstlps(thrtlrto, mach_to, engtype, engmode_toa, 'std', alt_to)
    instwf_torot = wf_torot(alt_to, mach_to, thstld, engtype, engmode_toa, thstlps_torot, instwf1)
    
    return (instwf_wu, instwf_toa, instwf_torot)

def wf_horzaccel(instwf_prev, alt, wngld, mach, thstlps, CLmax, velf, veli, engtype, eng_engmode, CDR=0, g0=9.81, gamma=1.4):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff rotation

    p 62
    """

    CL = 2*instwf_prev*(wngld/gamma/Pstd/get_atmos(alt, 'std')[1]/mach**2)
    CD0, K1, CD = get_drgplr(mach, CLmax, CL=CL)
    C1, C2 = get_engconst(engtype, eng_engmode)
    # print(CL, C1/mach+C2, CD)
    print(thstlps, instwf_prev)
    help1 = -(C1/mach+C2)*3600/spdsnd_std
    dvel = velf**2 - veli**2
    help3 = (CD+CDR)/CL/(instwf_prev/thstlps)/wngld
    # print(help1/3600*spdsnd_std)
    # print(dvel/2/g0/12/.0254)
    # print(help3, CD, CL, instwf_prev, thstlps)
    help2 = dvel/2/g0/(1-help3)
    instwf = np.exp(help1*help2)
    print(instwf)
    return instwf

def mission_anlys(thstld, wngld, thrtlrto, engtype, CLmax):
    """
    Notes
    -----
    Perform a weight analysis over the entire mission set
    """
    wngld = units.convert_pressure(wngld/144, 'Pa')
    alt_to = 2000
    instwfs_to = wf_takeoff(alt_to, thstld, wngld, thrtlrto, engtype, CLmax)
    instwf_to = instwfs_to[0]*instwfs_to[1]*instwfs_to[2]
    print(instwf_to)
    mach = .4410
    veli = 211.1*.0254*12
    velf = 812*.0254*12
    alt = 2000
    thstlps = get_thstlps(thrtlrto, mach, engtype, 'mil', 'trop', 2000)
    wf_horzaccel(instwf_to, alt, wngld, mach, thstlps, CLmax, velf, veli, engtype, 'mil')
    
    
def main():
    engtype = 'lbtf' # low bypass turbofan

    # constraint_anlys(engtype)

    # Aircraft intrinsic properties
    CLmax = 2
    # Chosen design point
    thstld = 1.25
    wngld = 64 # lbf/ft^2
    thrtlrto = 1.07

    # power_anlys(thstld, wngld, thrtlrto, engtype, CLmax)
    mission_anlys(thstld, wngld, thrtlrto, engtype, CLmax)

if __name__ == '__main__':
    # main()

    # thstlps = get_thstlps(1.07, .4410, 'lbtf', 'mil', 'std', 2000)
    # print(thstlps)
    hs = np.linspace(0, 30000, 4)/.0254/12
    t1 = get_atmos(hs, 'std')
    # t2 = get_atmos(30000, 'cold')
    # t3 = get_atmos(30000, 'hot')
    # t4 = get_atmos(30000, 'trop')

    # print(t1, t2, t3, t4)
    print(hs*.0254*12)
    print(t1)