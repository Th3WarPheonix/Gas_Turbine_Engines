
import numpy as np
from ambiance import Atmosphere as atmos
import Unit_Conversions as units
import matplotlib.pyplot as plt

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


def get_drag_polar(mach0:int, CLmax, vel_ratio=1.25, K2=0, sup=False):
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
    CD = CDmin + K'C^2 + K"(CL - CLmin)^2
    CD = K1CL^2 + K2CL + CDO
    K1 = K' + K"
    K2 = -2K"CLmin
    CD0 = CDmin + K"Cmin^2
    K' = induced drag
    K" = viscous drag
    Kprime = 1/pi/aspect_ratio/efficiency
    
    pg 37

    Returns
    ------
    0: CD0
    1: K1
    2: CD, drag coefficient

    Parameters
    ----------
    mach0 : freestream mach number
    CLmax : maximum coefficient of lift
    vel_ratio : velocity ratio V = vel_ratio*Vstall
    """

    if 0 <= mach0 <= 1:
        K1 = 0.18
    elif 1 < mach0 <= 2:
        K1 = (.18 - .36)/(1-2)*(mach0-1)+.18
    else:
        K1 = 0
        if not sup: print("Flight mach number not in range for turbine engines1")

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
        if not sup: print("Flight mach number not in range for turbine engines2")
        
    CL = CLmax/vel_ratio**2
    CD = K1*CL**2+K2*CL+CD0

    return (CD0, K1, CD)

def get_thrust_lapse(throttle_ratio:float, mach0, engine_type,
                     mode, alt:float=30000, gamma=1.4):
    """
    Notes
    -----
    Alpha is the installed full throttle thrust lapse, which depends on
    altitude, speed, and whether or not an afterburner is operating. i.e
    thrust = (alpha)(thrust_sea_level)

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

    if norm_Tt <= throttle_ratio:
        zeta = 0
    else:
        zeta = 1

    match engine_type:
        case 'hbtf': # high bypass turbofan
            thrust_lapse = norm_Pt*(1-0.49*np.sqrt(mach0) - 
                                    3*(norm_Tt-throttle_ratio)/
                                    (1.5+mach0)*zeta)
        case 'lbtf': # low bypass turbofan
            if mode == 'mil': # military power
                thrust_lapse = 0.6*norm_Pt*(1-3.8*(norm_Tt-throttle_ratio)/
                                            norm_Tt*zeta)
            elif mode == 'max': # maximum power
                thrust_lapse = norm_Pt*(1-3.5*(norm_Tt-throttle_ratio)/
                                        norm_Tt*zeta)
        case 'tj': # turbojet
            if mode == 'mil': # military power
                thrust_lapse = 0.8*norm_Pt*(1-.16*np.sqrt(mach0) - 
                                            25/norm_Tt*(norm_Tt-throttle_ratio)
                                            /(9+mach0)*zeta)
            elif mode == 'max': # maximum power
                thrust_lapse = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(mach0) -
                                         1.5*(norm_Tt-throttle_ratio)/
                                         norm_Tt*zeta)
        case 'tp': # turboprop
            thrust_lapse = norm_Pt # less than mach = 0.1
            thrust_lapse = norm_Pt*(1-0.96*(mach0-1)**0.25 - 3/8.13*(
                norm_Tt-throttle_ratio)/(mach0-0.1)*zeta)
    
    return thrust_lapse

def thrst_ld_takeoff(wing_load, densitys, CLmax, mach, alpha, beta, 
                        friction_coeff, vel_ratio_to, takeoff_distance, 
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
    CD = get_drag_polar(mach, CLmax, vel_ratio_to)[2]

    epsilon_to = CD + CDR + friction_coeff*CL
    third = beta/densitys/g0/epsilon_to
    second = rot_time*vel_ratio_to*np.sqrt(2*beta/densitys/CLmax)
    first = -(takeoff_distance-second*np.sqrt(wing_load))/third/wing_load
    thrust_load = (epsilon_to*vel_ratio_to**2/(1-np.exp(first))/CLmax+
                      friction_coeff)*beta/alpha

    return thrust_load

def thrst_ld_cruise(wing_load, alt, mach, CLmax, alpha, beta, 
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
    thrust_load = beta/alpha*(K1*beta/dyn_press*wing_load + K2 + 
                                 (CD0+CDR)/beta*dyn_press/wing_load) 

    return thrust_load

def thrst_ld_turn(wing_load, load_factor, alt, mach, CLmax, alpha, 
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
    thrust_load = beta/alpha*(first*wing_load + K2*load_factor + 
                                 second/wing_load) 

    return thrust_load

def thrst_ld_horz_accel(wing_load, deltamach, deltat, air_temp, 
                                 alt, mach, CLmax, alpha, beta, CDR=0, K2=0, 
                                 g0=9.81, gamma=1.4, R_gas=287.05):
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
    third = np.sqrt(gamma*R_gas*air_temp)*(deltamach)/g0/deltat
    thrust_load = beta/alpha*(first*wing_load + K2 + 
                                 second/wing_load + third) 

    return thrust_load

def thrst_ld_landing(wing_load, densitys, CLmax, CD, alpha, beta, 
                        friction_coeff, vel_ratio, distance, rotation_time=3, 
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
    epsilon = CD + CDR + friction_coeff*CL
    third = beta/densitys/g0/epsilon
    second = rotation_time*vel_ratio*np.sqrt(2*beta/densitys/CLmax)
    
    if alpha != 0:
        first = (distance-second*np.sqrt(wing_load))/third/wing_load
        thrust_load = (epsilon*vel_ratio**2/(np.exp(first)-1)/CLmax-
                          friction_coeff)*beta/alpha
    else:
        fourth = third*np.log(1+epsilon/friction_coeff/CLmax*vel_ratio**2)
        thrust_load = ((-second+np.sqrt(second**2+4*fourth*distance))/
                          2/fourth)**2

    return thrust_load

def excess_power(wing_load, thrust_load, load_factor, velocity, alt, 
                 mach, CLmax, alpha, beta, K2=0, CDR=0):
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
    CD0, K1, CD = get_drag_polar(mach, CLmax, sup=True)

    first  = alpha/beta*thrust_load
    second = -K1*load_factor**2*beta/dyn_press*wing_load
    third  = -K2*load_factor
    fourth = -(CD0+CDR)/beta/wing_load*dyn_press
    exs_pwr = velocity*(first+second+third+fourth)

    return exs_pwr 

def constraint_analysis(eng_type):
    """
    Notes
    -----
    Need to aid remining constraint equations form chapter 2 but
    currently not necessary for remainder of book walkthrough
    """
    
    wing_load1a = np.linspace(20, 120, 30)
    wing_load1 = units.convert_pressure(wing_load1a/144, 'Pa')
    
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

    alpha1 = get_thrust_lapse(TR, mach[0], eng_type, 'max', alt[0])
    alpha2 = get_thrust_lapse(TR, mach[1], eng_type, 'mil', alt[1])
    alpha3 = get_thrust_lapse(TR, mach[2], eng_type, 'max', alt[2])
    alpha4 = get_thrust_lapse(TR, mach[3], eng_type, 'max', alt[3])
    alpha5 = get_thrust_lapse(TR, mach[4], eng_type, 'max', alt[4])
    alpha7 = get_thrust_lapse(TR, mach[6], eng_type, 'max', alt[6])
    # All thrust loading equations have been validated
    thrst_ld1 = thrst_ld_takeoff(wing_load1, densitys, CLmax, mach[0], 
                                 alpha1, beta[0], muto, kto, takeoff_distance) 
    thrst_ld2 = thrst_ld_cruise(wing_load1, alt[1], mach[1], CLmax, alpha2,
                                beta[1]) 
    thrst_ld3 = thrst_ld_turn(wing_load1, load_turn1, alt[2], mach[2], 
                              CLmax, alpha3, beta[2]) 
    thrst_ld4 = thrst_ld_turn(wing_load1, load_turn2, alt[3], mach[3], 
                              CLmax, alpha4, beta[3]) 
    thrst_ld5 = thrst_ld_horz_accel(wing_load1, deltamach, deltat, air_temp,
                                    alt[4], mach[4], CLmax, alpha5, beta[4]) 
    thrst_ld6 = thrst_ld_landing(wing_load1, densitys, CLmax, landing_drag, 
                                 thrust_reverse, beta[5], mutd, ktd, 
                                 takeoff_distance)
    thrst_ld7 = thrst_ld_cruise(wing_load1, alt[6], mach[6], CLmax, alpha7, 
                                beta[6]) 
    
    thrst_ld6 = units.convert_pressure(thrst_ld6, "psi")*144

    plt.plot(wing_load1a, thrst_ld1, label="Takeoff", linestyle='--')
    plt.plot(wing_load1a, thrst_ld2, label="Cruise", linestyle='--')
    plt.plot(wing_load1a, thrst_ld3, label="Turn 1", linestyle='--')
    plt.plot(wing_load1a, thrst_ld4, label="Turn 2", linestyle='--')
    plt.plot(wing_load1a, thrst_ld5, label="Accel", linestyle='--')
    plt.vlines(thrst_ld6, .4, 1.6, label="Landing", linestyle='--')
    plt.plot(wing_load1a, thrst_ld7, label="Max Mach", linestyle='--')

    plt.legend()
    plt.ylim([.4, 1.6])
    plt.xlim([20, 120])
    # plt.show()
    plt.close()

def power_analysis(thrust_load, wing_load, thrtl_ratio, eng_type, CLmax):
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
    wing_load = units.convert_pressure(wing_load/144, 'Pa')
    N = 50
    velocity = np.linspace(100, 1900, N)*.0254*12 # ft/sec
    altitude = np.linspace(0, 60000, N) # ft
    power = np.zeros((N, N))
     
    for i, alt in enumerate(altitude):
        for j, vel in enumerate(velocity):
            mach = (vel/atmos(alt).speed_of_sound)[0]
            alpha = get_thrust_lapse(thrtl_ratio, mach, eng_type, 'mil', alt)
            power[i][j] = excess_power(wing_load, thrust_load, load_factor, 
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
    alpha = get_thrust_lapse(thrtl_ratio, mach, eng_type, 'mil', alt)
    pwr = excess_power(wing_load, thrust_load, load_factor, vel, alt, mach, 
                       CLmax, alpha, beta)
    print(pwr/.0254/12)

def empty_wieght_frac(takeoff_weight, type='fighter'):
    """
    Notes
    -----
    Preliminary results for empty wieght fraction, the ratio of empty
    aircrat wieght to takeoff weight (We/Wto).

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

def tsfc_initial(mach0, theta, eng_type, mode):
    """
    Notes
    -----
    Provide an inital estimate for the thrust specific fuel consumption
    units of 1/hour

    p 71
    """
    match eng_type:
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

def eng_const(eng_type, mode):
    """
    Notes
    -----
    Get engine type specific constants for mission analysis

    Parameters
    ----------
    eng_type:
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
    match eng_type:
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

def bcm_bca(beta, CLmax, wing_load, mach_crit=0.9, CDR=0, Pstd=101325, gamma=1.4):
    """
    Notes
    -----
    Obtain the best cruising altitude for subsonic speeds. The best
    cruise mach number is M=.9 after which divergence drag begins
    """
    CD0, K1, _ = get_drag_polar(mach_crit, CLmax)
    bca = 2*beta*wing_load/gamma/Pstd/mach_crit*2*np.sqrt((CD0+CDR)/K1)

    return bca

def wght_frac_warmup(alt, deltat, thrust_load, eng_type, mode, alpha, beta):
    """
    Notes Calculates wieght fraction of aircraft after warm-up

    p 67
    """
    theta = get_atmos(alt)[0]
    C1, _ = eng_const(eng_type, mode)

    wght_frac = 1 - C1/3600*np.sqrt(theta)*alpha/beta*thrust_load*deltat
    
    return wght_frac

def wght_frac_toaccel(alt, drag_coeff, fric_coeff, mach, mach_to, vel, thrust_load, wing_load, eng_type,
                      mode, alpha, beta, g0=9.81):
    """
    Notes Calculates wieght fraction of aircraft after takeoff
    acceleration

    p 63
    """
    theta = get_atmos(alt)[0]
    C1, C2 = eng_const(eng_type, mode)
    dyn_press = get_dyn_press(alt, mach_to)

    ttldrag_thrst = (drag_coeff*dyn_press/beta/wing_load+fric_coeff)*(
        beta/alpha/thrust_load)

    wght_frac = np.exp(-(C1+C2*mach)/3600*np.sqrt(theta)/g0*vel/
                        (1-ttldrag_thrst))

    return wght_frac

def wght_frac_torot(alt, mach_to, thrust_load, eng_type, mode, alpha, beta, 
                    rot_time=3):
    """
    Notes Calculates wieght fraction of aircraft after takeoff rotation

    p 68
    """
    theta = get_atmos(alt)[0]
    C1, C2 = eng_const(eng_type, mode)

    wght_frac = 1-(C1+C2*mach_to)/3600*np.sqrt(theta)*alpha/beta*(
        thrust_load*rot_time)

    return wght_frac

def get_takeoff_vel(beta, CLmax, densitys, wing_load, kto=1.2):
    """
    Notes
    -----
    Calculates the takeoff velocity
    
    p 29
    """
    vel_to = kto*np.sqrt(2*beta/densitys/CLmax*wing_load)
    return vel_to

def mission_analysis(thrust_load, wing_load, thrtl_ratio, eng_type, CLmax):
    """
    Notes
    -----
    Function to compare this codes results to book results
    """
    wing_load = units.convert_pressure(wing_load/144, 'Pa')
    alt_to = 2000

    mach_wu = 0
    beta_wu = 1
    deltat = 60 # s
    mode_wu = 'mil'
    alpha_wu = get_thrust_lapse(thrtl_ratio, mach_wu, eng_type, mode_wu, alt_to)
    beta_wu = wght_frac_warmup(alt_to, deltat, thrust_load, eng_type, mode_wu, alpha_wu, beta_wu)

    mach_to = .1819
    mach_toa = .1
    drag_coeff_toa = .3612
    fric_coeff_toa = .05
    mode_toa = 'max'
    densitys = get_atmos(alt_to)[4]
    vel_to = get_takeoff_vel(beta_wu, CLmax, densitys, wing_load)
    alpha_toa = get_thrust_lapse(thrtl_ratio, mach_toa, eng_type, mode_toa, alt_to)
    beta_toa = wght_frac_toaccel(alt_to, drag_coeff_toa, fric_coeff_toa, mach_toa, mach_to, vel_to, thrust_load, wing_load, eng_type,
                      mode_toa, alpha_toa, beta_wu)
    
    beta1 = beta_toa*beta_wu
    alpha_torot = get_thrust_lapse(thrtl_ratio, mach_to, eng_type, mode_toa, alt_to)
    beta_torot = wght_frac_torot(alt_to, mach_to, thrust_load, eng_type, mode_toa, alpha_torot, beta1)

    beta2 = beta1*beta_torot
    
    print(beta2)
    
    
def main():
    eng_type = 'lbtf' # low bypass turbofan

    # constraint_analysis(eng_type)

    # Aircraft intrinsic properties
    CLmax = 2
    # Chosen design point
    thrust_load = 1.25
    wing_load = 64 # lbf/ft^2
    throttle_ratio = 1.07

    # power_analysis(thrust_load, wing_load, throttle_ratio, eng_type, CLmax)
    mission_analysis(thrust_load, wing_load, throttle_ratio, eng_type, CLmax)

if __name__ == '__main__':
    main()