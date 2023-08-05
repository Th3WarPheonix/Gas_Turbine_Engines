
import numpy as np
from ambiance import Atmosphere as atmos
import Unit_Conversions as units
import matplotlib.pyplot as plt

Pstd = 101325 # Pa
Tstd = 288.15 # K
densitystd = 1.225 # kg/m^3
sound_speed_std = 340.3 # m/s

"""
Variables 
---------
variable name : actual name : quantity/description : in book

wngld : wing loading : takeoff_weight/wing_planform_area : Wto/S
thstld : thrust laoding : thrust/takeoff_weight : T/Wto
engtype : engine type : tj, lbtf, tp, etc : 
thrtlrto : throttle ratio :  : TR
thstlps : installed thrust lapse : thrust/sea_level_thrust: alpha
instwf : instantaneous weight fraction : weight/takeoff_weight : beta

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

Abbreviations and tags
----------------------
BCA : best cruise altitude
BCM : best cruise altitude
SLS : sea level static
dry : afterburner off
wet : afterburner on
anlys : analysis
horz : horizontal
accel : acceleration
ln : landing
to : takeoff
td : touchdown
obs : obstacle

"""

def mass_flow_parameter(mach, gamma, Rgas):
    """
    Returns
    -------
    mfp : mass flow parameter m*sqrt(Tt)/Pt/A
    """
    exponent = (gamma+1)/2/(1-gamma)
    mfp = mach*np.sqrt(gamma/Rgas)*((1+(gamma-1)/2)*mach**2)**exponent
    return mfp

def get_atmos(altitude):
    """
    Notes
    -----
    Altitude in feet. Can be a numpy array
    
    Returns
    -------
    0: theta Tt/Tstd
    1: delta Pt/Pstd
    2: temperature
    3: pressure
    4: density
    """
    temp = atmos(altitude*.0254*12).temperature[0]
    press = atmos(altitude*.0254*12).pressure[0]
    density = atmos(altitude*.0254*12).pressure[0]
    
    theta = temp/288.15
    delta = press/101325
    if altitude == 2000:
        theta = 1.0796
        delta = 0.9298
        density = delta/theta*1.225
    return (theta, delta, temp, press, density)

def get_dyn_press(alt, mach, gamma=1.4):
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

    delta= get_atmos(alt)[1]
    return delta*101325*mach**2*gamma/2


def get_drag_polar(mach0:int, CLmax, velrt=1.25, K2=0, CL=0, suppress=False):
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

def get_thstlps(thrtlrto:float, mach0, engine_type,
                     mode, alt:float=30000, gamma=1.4):
    """
    Notes
    -----
    Alpha is the installed full throttle thrust lapse, which depends on
    altitude, speed, and whether or not an afterburner is operating. i.e
    alpha = thrust/thrust_sea_level

    pg 38
    
    Parameters
    ----------
    engine_type : 
        hbtf : High bypass turbofan
        lbtf : Low bypass turbofan
        tj : Turbojet
        tp : Turboprop
    mode : 
        mil : military power
        max : max power
    alt : altitude in feet
    """
    
    if alt == 2000:
        theta = 1.0796
        delta = .9298
    else:
        theta = get_atmos(alt)[0]
        delta = get_atmos(alt)[1]

    # normalized total temperature and total pressure
    norm_Tt = theta*(1+(gamma-1)/2*mach0**2)
    norm_Pt = delta*(1+(gamma-1)/2*mach0**2)**(gamma/(gamma-1))

    if norm_Tt <= thrtlrto:
        zeta = 0
    else:
        zeta = 1

    match engine_type:
        case 'hbtf': # high bypass turbofan
            thstlps = norm_Pt*(1-0.49*np.sqrt(mach0) - 
                                    3*(norm_Tt-thrtlrto)/
                                    (1.5+mach0)*zeta)
        case 'lbtf': # low bypass turbofan
            if mode == 'mil': # military power
                thstlps = 0.6*norm_Pt*(1-3.8*(norm_Tt-thrtlrto)/
                                            norm_Tt*zeta)
            elif mode == 'max': # maximum power
                thstlps = norm_Pt*(1-3.5*(norm_Tt-thrtlrto)/
                                        norm_Tt*zeta)
        case 'tj': # turbojet
            if mode == 'mil': # military power
                thstlps = 0.8*norm_Pt*(1-.16*np.sqrt(mach0) - 
                                            25/norm_Tt*(norm_Tt-thrtlrto)
                                            /(9+mach0)*zeta)
            elif mode == 'max': # maximum power
                thstlps = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(mach0) -
                                         1.5*(norm_Tt-thrtlrto)/
                                         norm_Tt*zeta)
        case 'tp': # turboprop
            thstlps = norm_Pt # less than mach = 0.1
            thstlps = norm_Pt*(1-0.96*(mach0-1)**0.25 - 3/8.13*(
                norm_Tt-thrtlrto)/(mach0-0.1)*zeta)
    
    return thstlps


def get_takeoff_vel(beta, CLmax, densitys, wngld, kto=1.2):
    """
    Notes
    -----
    Calculates the takeoff velocity
    
    p 29
    """
    vel_to = kto*np.sqrt(2*beta/densitys/CLmax*wngld)
    return vel_to


def thstld_takeoff(wngld, densitys, CLmax, mach, alpha, beta, 
                        frccf_ln, velrt_to, takeoff_distance, 
                        rot_time=3, CL=0, CDR=0, g0=9.81):
    """
    Notes
    -----
    Thrust loading during takeoff acceleration up to rotation
    instantaneous weight is
    given by
    W = beta*WTO
    beta, instantaneous weight fraction, depends on how much fuel has
    been consumed and payload delivered; more AED2e p24

    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    standard rotation time is 3 secs
    epsilon is the sum of the drags on the aircraft
    kto is the ratio of takeoff veloicty to stall velocity

    pg 34
    """
    CD = get_drag_polar(mach, CLmax, velrt_to)[2]

    epsilon_to = CD + CDR + frccf_ln*CL
    third = beta/densitys/g0/epsilon_to
    second = rot_time*velrt_to*np.sqrt(2*beta/densitys/CLmax)
    first = -(takeoff_distance-second*np.sqrt(wngld))/third/wngld
    thstld = (epsilon_to*velrt_to**2/(1-np.exp(first))/CLmax+
                      frccf_ln)*beta/alpha

    return thstld

def thstld_cruise(wngld, alt, mach, CLmax, alpha, beta, 
                       CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during cruise
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 26
    """
    dyn_press = get_dyn_press(alt, mach)
    CD0, K1, CD = get_drag_polar(mach, CLmax)
    thstld = beta/alpha*(K1*beta/dyn_press*wngld + K2 + 
                                 (CD0+CDR)/beta*dyn_press/wngld) 

    return thstld

def thstld_turn(wngld, load_factor, alt, mach, CLmax, alpha, 
                     beta, CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 26
    """
    dyn_press = get_dyn_press(alt, mach)
    CD0, K1, CD = get_drag_polar(mach, CLmax)

    first  = K1*load_factor**2*beta/dyn_press
    second = (CD0+CDR)/beta*dyn_press
    thstld = beta/alpha*(first*wngld + K2*load_factor + 
                                 second/wngld) 

    return thstld

def thstld_horzaccel(wngld, deltamach, deltat, air_temp, alt, mach, CLmax, 
                     alpha, beta, CDR=0, K2=0, g0=9.81, gamma=1.4, Rgas=287.05):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S

    pg 27
    """
    dyn_press = get_dyn_press(alt, mach)
    CD0, K1, CD = get_drag_polar(mach, CLmax)

    first  = K1*beta/dyn_press
    second = (CD0+CDR)/beta*dyn_press
    third = np.sqrt(gamma*Rgas*air_temp)*(deltamach)/g0/deltat
    thstld = beta/alpha*(first*wngld + K2 + 
                                 second/wngld + third) 

    return thstld

def thstld_landing(wngld, densitys, CLmax, CD, alpha, beta, frccf_ln, velrt, distance, rotation_time=3, 
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
    third = beta/densitys/g0/epsilon
    second = rotation_time*velrt*np.sqrt(2*beta/densitys/CLmax)
    
    if alpha != 0:
        first = (distance-second*np.sqrt(wngld))/third/wngld
        thstld = (epsilon*velrt**2/(np.exp(first)-1)/CLmax-
                          frccf_ln)*beta/alpha
    else:
        fourth = third*np.log(1+epsilon/frccf_ln/CLmax*velrt**2)
        thstld = ((-second+np.sqrt(second**2+4*fourth*distance))/
                          2/fourth)**2

    return thstld

def excess_power(wngld, thstld, load_factor, velocity, alt, mach, CLmax, alpha,
                 beta, K2=0, CDR=0):
    """
    Notes
    -----
    Calculates excess power for the aircraft at the defined thrust
    loading, wing loading, throttle ratio for any beta, load factor
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
    dyn_press = get_dyn_press(alt, mach)
    CD0, K1, CD = get_drag_polar(mach, CLmax, suppress=True)

    first  = alpha/beta*thstld
    second = -K1*load_factor**2*beta/dyn_press*wngld
    third  = -K2*load_factor
    fourth = -(CD0+CDR)/beta/wngld*dyn_press
    exs_pwr = velocity*(first+second+third+fourth)

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
    beta = [1, .78, .78, .78, .78, .56, .78]
    takeoff_distance = 1500*12*.0254 # m
    mach = [.1, 1.5, 1.6, .9, 1.2, 0, 1.8]
    alt = [2000, 30000, 30000, 30000, 30000, 0, 40000]

    load_turn1 = load_turn2 = 5 # g

    deltamach = .8
    deltat = 50

    air_temp = get_atmos(30000)[0]*288.15
    kto = 1.2
    ktd = 1.15
    muto = .05
    mutd = .18
    CLmax = 2

    landing_drag = .8123
    thrust_reverse = 0

    alpha1 = get_thstlps(TR, mach[0], engtype, 'max', alt[0])
    alpha2 = get_thstlps(TR, mach[1], engtype, 'mil', alt[1])
    alpha3 = get_thstlps(TR, mach[2], engtype, 'max', alt[2])
    alpha4 = get_thstlps(TR, mach[3], engtype, 'max', alt[3])
    alpha5 = get_thstlps(TR, mach[4], engtype, 'max', alt[4])
    alpha7 = get_thstlps(TR, mach[6], engtype, 'max', alt[6])
    # All thrust loading equations have been validated
    thstld1 = thstld_takeoff(wngld1, densitys, CLmax, mach[0], alpha1, 
                                 beta[0], muto, kto, takeoff_distance) 
    thstld2 = thstld_cruise(wngld1, alt[1], mach[1], CLmax, alpha2,
                                beta[1]) 
    thstld3 = thstld_turn(wngld1, load_turn1, alt[2], mach[2], 
                              CLmax, alpha3, beta[2]) 
    thstld4 = thstld_turn(wngld1, load_turn2, alt[3], mach[3], CLmax, 
                              alpha4, beta[3]) 
    thstld5 = thstld_horzaccel(wngld1, deltamach, deltat, air_temp,
                                    alt[4], mach[4], CLmax, alpha5, beta[4]) 
    thstld6 = thstld_landing(wngld1, densitys, CLmax, landing_drag, 
                                 thrust_reverse, beta[5], mutd, ktd, 
                                 takeoff_distance)
    thstld7 = thstld_cruise(wngld1, alt[6], mach[6], CLmax, alpha7, 
                                beta[6]) 
    
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
    # plt.show()
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
    beta = .97
    wngld = units.convert_pressure(wngld/144, 'Pa')
    N = 50
    velocity = np.linspace(100, 1900, N)*.0254*12 # ft/sec
    altitude = np.linspace(0, 60000, N) # ft
    power = np.zeros((N, N))
     
    for i, alt in enumerate(altitude):
        for j, vel in enumerate(velocity):
            mach = (vel/atmos(alt).speed_of_sound)[0]
            alpha = get_thstlps(thrtlrto, mach, engtype, 'mil', alt)
            power[i][j] = excess_power(wngld, thstld, load_factor, 
                                 vel, alt, mach, CLmax, alpha, beta)
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
    alpha = get_thstlps(thrtlrto, mach, engtype, 'mil', alt)
    pwr = excess_power(wngld, thstld, load_factor, vel, alt, mach, 
                       CLmax, alpha, beta)
    print(pwr/.0254/12)

def empty_weight_frac(takeoff_weight, type='fighter'):
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
            we_frac = 2.34*takeoff_weight**-0.13
        case 'cargo':
            we_frac = 1.26*takeoff_weight**-0.08
        case 'passenger':
            we_frac = 1.02*takeoff_weight**-0.06
        case 'twin tprop':
            we_frac = 0.96*takeoff_weight**-0.05

    return we_frac

def tsfc_initial(mach0, theta, engtype, mode):
    """
    Notes
    -----
    Provide an inital estimate for the thrust specific fuel consumption
    units of 1/hour

    p 71
    """
    match engtype:
        case 'hbtf': # high bypass turbofan
            tsfc = (0.45+0.54*mach0)*np.sqrt(theta)
        case 'lbtf': # low bypass turbofan
            if mode == 'mil': # military power
                tsfc = (0.9+0.30*mach0)*np.sqrt(theta)
            elif mode == 'max': # maximum power
                tsfc = (1.6+0.27*mach0)*np.sqrt(theta)
        case 'tj': # turbojet
            if mode == 'mil': # military power
                tsfc = (1.1+0.30*mach0)*np.sqrt(theta)
            elif mode == 'max': # maximum power
                tsfc = (1.5+0.23*mach0)*np.sqrt(theta)
        case 'tp': # turboprop
            tsfc = (0.18+0.8*mach0)*np.sqrt(theta)

    return tsfc

def get_engconst(engtype, mode):
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
    mode:
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
            if mode == 'mil': # military power
                consts = [0.9, 0.30]
            elif mode == 'max': # maximum power
                consts = [1.6, 0.27]
        case 'tj': # turbojet
            if mode == 'mil': # military power
                consts = [1.1, 0.30]
            elif mode == 'max': # maximum power
                consts = [1.5, 0.23]
        case 'tp': # turboprop
            consts = [0.18, 0.8]
    
    return consts

def get_bcmbca(beta, CLmax, wngld, mach_crit=0.9, CDR=0, Pstd=101325, gamma=1.4):
    """
    Notes
    -----
    Obtain the best cruising altitude for subsonic speeds. The best
    cruise mach number is M=.9 after which divergence drag begins
    """
    CD0, K1, _ = get_drag_polar(mach_crit, CLmax)
    bca = 2*beta*wngld/gamma/Pstd/mach_crit*2*np.sqrt((CD0+CDR)/K1)

    return bca

def wf_warmup(alt, deltat, thstld, engtype, mode, alpha, beta):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after warm-up

    p 67
    """
    theta = get_atmos(alt)[0]
    C1, _ = get_engconst(engtype, mode)

    wght_frac = 1 - C1/3600*np.sqrt(theta)*alpha/beta*thstld*deltat
    
    return wght_frac

def wf_toaccel(alt, drag_coeff, fric_coeff, mach, mach_to, vel, thstld, wngld, engtype,
                      mode, alpha, beta, g0=9.81):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff
    acceleration

    p 63
    """
    theta = get_atmos(alt)[0]
    C1, C2 = get_engconst(engtype, mode)
    dyn_press = get_dyn_press(alt, mach_to)

    ttldrag_thrst = (drag_coeff*dyn_press/beta/wngld+fric_coeff)*(
        beta/alpha/thstld)

    wght_frac = np.exp(-(C1+C2*mach)/3600*np.sqrt(theta)/g0*vel/
                        (1-ttldrag_thrst))

    return wght_frac

def wf_torot(alt, mach_to, thstld, engtype, mode, alpha, beta, 
                    rot_time=3):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff rotation

    p 68
    """
    theta = get_atmos(alt)[0]
    C1, C2 = get_engconst(engtype, mode)

    wght_frac = 1-(C1+C2*mach_to)/3600*np.sqrt(theta)*alpha/beta*(
        thstld*rot_time)

    return wght_frac

def wf_takeoff(alt_to, thstld, wngld, thrtlrto, engtype, CLmax):
    """
    Notes
    -----
    Function to compare this codes results to book results
    """
    
    mach_wu = 0
    beta_wu = 1
    deltat = 60 # s
    mode_wu = 'mil'
    alpha_wu = get_thstlps(thrtlrto, mach_wu, engtype, mode_wu, alt_to)
    beta_wu = wf_warmup(alt_to, deltat, thstld, engtype, mode_wu, alpha_wu, beta_wu)

    mach_to = .1819
    mach_toa = .1
    drag_coeff_toa = .3612
    fric_coeff_toa = .05
    mode_toa = 'max'
    densitys = get_atmos(alt_to)[4]
    vel_to = get_takeoff_vel(beta_wu, CLmax, densitys, wngld)
    alpha_toa = get_thstlps(thrtlrto, mach_toa, engtype, mode_toa, alt_to)
    beta_toa = wf_toaccel(alt_to, drag_coeff_toa, fric_coeff_toa, mach_toa, mach_to, vel_to, thstld, wngld, engtype,
                      mode_toa, alpha_toa, beta_wu)
    
    beta1 = beta_toa*beta_wu
    alpha_torot = get_thstlps(thrtlrto, mach_to, engtype, mode_toa, alt_to)
    beta_torot = wf_torot(alt_to, mach_to, thstld, engtype, mode_toa, alpha_torot, beta1)
    
    return (beta_wu, beta_toa, beta_torot)


def wf_horzaccel(beta_prev, alt, wngld, mach, alpha, CLmax, velf, veli, engtype, eng_mode, CDR=0, g0=9.81, gamma=1.4):
    """
    Notes
    -----
    Calculates weight fraction of aircraft after takeoff rotation

    p 62
    """

    CL = 2*beta_prev*(wngld/gamma/Pstd/get_atmos(alt)[1]/mach**2)
    CD0, K1, CD = get_drag_polar(mach, CLmax, CL=CL)
    C1, C2 = get_engconst(engtype, eng_mode)

    first = -(C1/mach+C2)*3600/sound_speed_std
    dvel = velf**2 - veli**2
    third = (CD+CDR)/CL/(beta_prev/alpha)/wngld
    print(first/3600*sound_speed_std)
    print(dvel/2/g0/12/.0254)
    print(third, CD, CL, beta_prev, alpha)
    second = dvel/2/g0/(1-third)
    beta = np.exp(first*second)
    print(beta)
    return beta

def mission_anlys(thstld, wngld, thrtlrto, engtype, CLmax):
    """
    Notes
    -----
    Perform a weight anlys over the entire mission set
    """
    wngld = units.convert_pressure(wngld/144, 'Pa')
    alt_to = 2000
    betas_to = wf_takeoff(alt_to, thstld, wngld, thrtlrto, engtype, CLmax)
    beta_to = betas_to[0]*betas_to[1]*betas_to[2]

    mach = .4410
    veli = 211.1*.0254*12
    velf = 812*.0254*12
    alt = 2000
    alpha = get_thstlps(thrtlrto, mach, engtype, 'mil')
    wf_horzaccel(beta_to, alt, wngld, mach, alpha, CLmax, velf, veli, engtype, 'mil')
    
    
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
    main()