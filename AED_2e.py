
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
    2: 
    """
    theta = atmos(altitude*.0254*12).temperature[0]/288.15
    delta = atmos(altitude*.0254*12).pressure[0]/101325

    return (theta, delta)

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

    theta, delta = get_atmos(alt)
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

def get_thrust_lapse(throttle_ratio:float, mach0, engine_type:float,
                     alt:float=30000, mode:int=0, gamma=1.4):
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
        0 : High Bypass Turbofan
        1 : Low Bypass Turbofan, mixed flow
        2 : Turbojet 
        3 : Turboprop
    mode : 
        0 : military power
        1 : max power
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
        case 0: # high bypass turbofan
            thrust_lapse = norm_Pt*(1-0.49*np.sqrt(mach0) - 
                                    3*(norm_Tt-throttle_ratio)/
                                    (1.5+mach0)*zeta)
        case 1: # low bypass turbofan
            if mode == 0: # military power
                thrust_lapse = 0.6*norm_Pt*(1-3.8*(norm_Tt-throttle_ratio)/
                                            norm_Tt*zeta)
            elif mode == 1: # maximum power
                thrust_lapse = norm_Pt*(1-3.5*(norm_Tt-throttle_ratio)/
                                        norm_Tt*zeta)
        case 2: # turbojet
            if mode == 0: # military power
                thrust_lapse = 0.8*norm_Pt*(1-.16*np.sqrt(mach0) - 
                                            25/norm_Tt*(norm_Tt-throttle_ratio)
                                            /(9+mach0)*zeta)
            elif mode == 1: # maximum power
                thrust_lapse = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(mach0) -
                                         1.5*(norm_Tt-throttle_ratio)/
                                         norm_Tt*zeta)
        case 3: # turboprop
            thrust_lapse = norm_Pt # less than mach = 0.1
            thrust_lapse = norm_Pt*(1-0.96*(mach0-1)**0.25 - 3/8.13*(
                norm_Tt-throttle_ratio)/(mach0-0.1)*zeta)
    
    return thrust_lapse

def thrst_ld_takeoff(wing_loading, densitys, CLmax, mach, alpha, beta, 
                        friction_coeff, vel_ratio_to, takeoff_distance, 
                        rotation_time=3, CL=0, CDR=0, g0=9.81):
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
    second = rotation_time*vel_ratio_to*np.sqrt(2*beta/densitys/CLmax)
    first = -(takeoff_distance-second*np.sqrt(wing_loading))/third/wing_loading
    thrust_loading = (epsilon_to*vel_ratio_to**2/(1-np.exp(first))/CLmax+
                      friction_coeff)*beta/alpha

    return thrust_loading

def thrst_ld_cruise(wing_loading, alt, mach, CLmax, alpha, beta, 
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
    thrust_loading = beta/alpha*(K1*beta/dyn_press*wing_loading + K2 + 
                                 (CD0+CDR)/beta*dyn_press/wing_loading) 

    return thrust_loading

def thrst_ld_turn(wing_loading, load_factor, alt, mach, CLmax, alpha, 
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
    thrust_loading = beta/alpha*(first*wing_loading + K2*load_factor + 
                                 second/wing_loading) 

    return thrust_loading

def thrst_ld_horz_accel(wing_loading, deltamach, deltat, air_temp, 
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
    thrust_loading = beta/alpha*(first*wing_loading + K2 + 
                                 second/wing_loading + third) 

    return thrust_loading

def thrst_ld_landing(wing_loading, densitys, CLmax, CD, alpha, beta, 
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
        first = (distance-second*np.sqrt(wing_loading))/third/wing_loading
        thrust_loading = (epsilon*vel_ratio**2/(np.exp(first)-1)/CLmax-
                          friction_coeff)*beta/alpha
    else:
        fourth = third*np.log(1+epsilon/friction_coeff/CLmax*vel_ratio**2)
        thrust_loading = ((-second+np.sqrt(second**2+4*fourth*distance))/
                          2/fourth)**2

    return thrust_loading

def excess_power(wing_loading, thrust_loading, load_factor, velocity, alt, 
                 mach, CLmax, alpha, beta, K2=0, CDR=0):
    """
    Notes
    -----
    NOT VALIDATED
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

    first  = alpha/beta*thrust_loading
    second = -K1*load_factor**2*beta/dyn_press*wing_loading
    third  = -K2*load_factor
    fourth = -(CD0+CDR)/beta/wing_loading*dyn_press
    exs_pwr = velocity*(first+second+third+fourth)

    return exs_pwr 

def constraint_analysis():
    
    wing_loading1a = np.linspace(20, 120, 30)
    wing_loading1 = units.convert_pressure(wing_loading1a/144, 'Pa')
    
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

    alpha1 = get_thrust_lapse(TR, mach[0], 1, alt[0], 1)
    alpha2 = get_thrust_lapse(TR, mach[1], 1, alt[1], 0)
    alpha3 = get_thrust_lapse(TR, mach[2], 1, alt[2], 1)
    alpha4 = get_thrust_lapse(TR, mach[3], 1, alt[3], 1)
    alpha5 = get_thrust_lapse(TR, mach[4], 1, alt[4], 1)
    alpha7 = get_thrust_lapse(TR, mach[6], 1, alt[6], 1)
    # All thrust loading equations have been validated
    thrst_ld1 = thrst_ld_takeoff(wing_loading1, densitys, CLmax, mach[0], 
                                 alpha1, beta[0], muto, kto, takeoff_distance) 
    thrst_ld2 = thrst_ld_cruise(wing_loading1, alt[1], mach[1], CLmax, alpha2,
                                beta[1]) 
    thrst_ld3 = thrst_ld_turn(wing_loading1, load_turn1, alt[2], mach[2], 
                              CLmax, alpha3, beta[2]) 
    thrst_ld4 = thrst_ld_turn(wing_loading1, load_turn2, alt[3], mach[3], 
                              CLmax, alpha4, beta[3]) 
    thrst_ld5 = thrst_ld_horz_accel(wing_loading1, deltamach, deltat, air_temp,
                                    alt[4], mach[4], CLmax, alpha5, beta[4]) 
    thrst_ld6 = thrst_ld_landing(wing_loading1, densitys, CLmax, landing_drag, 
                                 thrust_reverse, beta[5], mutd, ktd, 
                                 takeoff_distance)
    thrst_ld7 = thrst_ld_cruise(wing_loading1, alt[6], mach[6], CLmax, alpha7, 
                                beta[6]) 
    
    thrst_ld6 = units.convert_pressure(thrst_ld6, "psi")*144

    plt.plot(wing_loading1a, thrst_ld1, label="Takeoff", linestyle='--')
    plt.plot(wing_loading1a, thrst_ld2, label="Cruise", linestyle='--')
    plt.plot(wing_loading1a, thrst_ld3, label="Turn 1", linestyle='--')
    plt.plot(wing_loading1a, thrst_ld4, label="Turn 2", linestyle='--')
    plt.plot(wing_loading1a, thrst_ld5, label="Accel", linestyle='--')
    plt.vlines(thrst_ld6, .4, 1.6, label="Landing", linestyle='--')
    plt.plot(wing_loading1a, thrst_ld7, label="Max Mach", linestyle='--')

    plt.legend()
    plt.ylim([.4, 1.6])
    plt.xlim([20, 120])
    # plt.show()
    plt.close()

def power_analysis():
    # Chosen design point
    thrust_loading = 1.25
    wing_loading = 64 # lbf/ft^2
    throttle_ratio = 1.07

    load_factor = 1
    beta = .97
    CLmax = 2
    wing_loading = units.convert_pressure(wing_loading/144, 'Pa')
    N = 100
    velocity = np.linspace(100, 2000, N)*.0254*12 # ft/sec
    altitude = np.linspace(0, 70000, N) # ft
    power = np.zeros((N, N))
    mach = velocity/atmos(altitude).speed_of_sound
   
    for i, alt in enumerate(altitude):
        alpha = get_thrust_lapse(throttle_ratio, mach[i], 1, alt, 0)
        for j, vel in enumerate(velocity):
            power[i][j] = excess_power(wing_loading, thrust_loading, load_factor, 
                                 vel, alt, mach[i], CLmax, alpha, beta)
    power = power/12/.0254
    
    levels = np.arange(0, 600+1, 50)
    x = np.linspace(13, 5, 10)
    y = np.linspace(-7, 5, 10)

    zz2 = np.zeros((10, 10))
    for i, xs in enumerate(x):
        for j, ys in enumerate(y):
            zz2[j][i] = np.sqrt(xs**2 + ys**2)


    plt.contourf(velocity/12/.0254, altitude, power, levels)
    plt.colorbar()
    plt.savefig('pic.png')

    alt = 36000
    vel = 1600*12*.0254
    mach = vel/(atmos(alt).speed_of_sound)[0]
    alpha = get_thrust_lapse(throttle_ratio, mach, 1, alt, 0)
    pwr = excess_power(wing_loading, thrust_loading, load_factor, vel, alt, mach, CLmax, alpha, beta)
    print(pwr/.0254/12)

if __name__ == '__main__':
    # constraint_analysis()
    power_analysis()
