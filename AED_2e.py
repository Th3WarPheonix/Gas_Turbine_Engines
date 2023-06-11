
import numpy as np
from ambiance import Atmosphere as atmos
import Unit_Conversions as units
import matplotlib.pyplot as plt

def get_atmos(altitude):
    """Altitude in feet
    
    Returns
    -------
    0: theta Tt/Tstd
    1: delta Pt/Pstd
    """
    theta = atmos(altitude*.0254*12).temperature[0]/288.15
    delta = atmos(altitude*.0254*12).pressure[0]/101325

    return (theta, delta)

def get_dyn_press(alt, mach, gamma=1.4):
    """
    Returns
    -------
    dynamic pressure in Pascals

    Parameters
    ----------
    alt : altitude in feet
    """

    theta, delta = get_atmos(alt)
    return delta*101325*mach**2*gamma/2


def drag_polar(mach0:int, CLmax, kto, K2=0):
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

    Parameters
    ----------
    
    """

    if 0 <= mach0 <= 1:
        K1 = 0.18
    elif 1 < mach0 <= 2:
        K1 = (.18 - .36)/(1-2)*(mach0-1)+.18
    else:
        print("Flight mach number not in range for turbine engines")

    if 0 <= mach0 <= 0.6:
        CD0 = .014
    elif  0.6 < mach0 <= 0.8:
        CD0 = (.014 - .016)/(.6-.8)*(mach0-0.6)+.014
    elif 0.8 < mach0 <= 1.2:
        CD0 = (.016 - .026)/(.8-1.2)*(mach0-0.8)+.016
    elif 1.2 < mach0 <= 1.4:
        CD0 = (.026 - .028)/(1.2-1.4)*(mach0-1.2)+.026
    elif 1.4 < mach0 <= 2:
        CD0 = .028
        
    CL = CLmax/kto**2
    CD = K1*CL**2+K2*CL+CD0

    return CD, CD0, K1

def engine_thrust_lapse(throttle_ratio, mach0, engine_type, alt=30, mode=0, gamma=1.4):
    """
    Notes
    -----
    alpha is the installed full throttle thrust lapse, which depends on altitude, speed,
    and whether or not an afterburner is operating
    Thrust = alpha*thrustsl
    
    Parameters
    ----------
    engine_type : 0 HBTF : 1 LBTF-MF : 2 TJ : 3 TP
    mode : 0 military power : 1 max power
    alt : in thousands of feet
    """
    # TURN THIS INTO A LOOK UP TABLE
    # normalized total temperature and total pressure

    if alt == 2:
        theta = 1.0796
        delta = .9298
    else:
        theta = get_atmos(alt*1000)[0]
        delta = get_atmos(alt*1000)[1]

    norm_Tt = theta*(1+(gamma-1)/2*mach0**2)
    norm_Pt = delta*(1+(gamma-1)/2*mach0**2)**(gamma/(gamma-1))

    if norm_Tt <= throttle_ratio:
        zeta = 0
    else:
        zeta = 1

    match engine_type:
        case 0: # high bypass turbofan
            thrust_lapse = norm_Pt*(1-0.49*np.sqrt(mach0) - 3*(norm_Tt-throttle_ratio)/(1.5+mach0)*zeta)
        case 1: # low bypass turbofan
            if mode == 0: # military power
                thrust_lapse = 0.6*norm_Pt*(1-3.8*(norm_Tt-throttle_ratio)/norm_Tt*zeta)
            elif mode == 1: # maximum power
                thrust_lapse = norm_Pt*(1-3.5*(norm_Tt-throttle_ratio)/norm_Tt*zeta)
        case 2: # turbojet
            if mode == 0: # military power
                thrust_lapse = 0.8*norm_Pt*(1-.16*np.sqrt(mach0) - 25/norm_Tt*(norm_Tt-throttle_ratio)/(9+mach0)*zeta)
            elif mode == 1: # maximum power
                thrust_lapse = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(mach0) - 1.5*(norm_Tt-throttle_ratio)/norm_Tt*zeta)
        case 3: # turboprop
            thrust_lapse = norm_Pt # less than mach = 0.1
            thrust_lapse = norm_Pt*(1-0.96*(mach0-1)**0.25 - 3/8.13*(norm_Tt-throttle_ratio)/(mach0-0.1)*zeta)

    # print(thrust_lapse, theta, norm_Tt, delta, norm_Pt)
    
    return thrust_lapse

def thrust_loading_takeoff(wing_loading, densitys, CLmax, CD, alpha, beta, friction_coeff, vel_ratio_to, takeoff_distance, rotation_time=3, CL=0, CDR=0, g0=9.81):
    """
    Notes
    -----
    Thrust loading during takeoff acceleration up to rotation
    instantaneous weight is
    given by
    W = beta*WTO
    beta, instantaneous weight fraction, depends on how much fuel has been consumed and payload delivered; more AED2e p24

    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    standard rotation time is 3 secs
    epsilon is the sum of the drags on the aircraft
    kto is the ratio of takeoff veloicty to stall velocity
    """

    epsilon_to = CD + CDR + friction_coeff*CL
    third = beta/densitys/g0/epsilon_to
    second = rotation_time*vel_ratio_to*np.sqrt(2*beta/densitys/CLmax)
    first = -(takeoff_distance-second*np.sqrt(wing_loading))/third/wing_loading
    thrust_loading = (epsilon_to*vel_ratio_to**2/(1-np.exp(first))/CLmax+friction_coeff)*beta/alpha

    return thrust_loading

def thrust_loading_cruise(wing_loading, dyn_press, CD0, K1, alpha, beta, CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during cruise
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    """
    thrust_loading = beta/alpha*(K1*beta/dyn_press*wing_loading + K2 + (CD0+CDR)/beta*dyn_press/wing_loading) 

    return thrust_loading

def thrust_loading_turn(wing_loading, load_factor, dyn_press, CD0, K1, alpha, beta, CDR=0, K2=0):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    """

    first  = K1*load_factor**2*beta/dyn_press
    second = (CD0+CDR)/beta*dyn_press
    thrust_loading = beta/alpha*(first*wing_loading + K2*load_factor + second/wing_loading) 

    return thrust_loading

def thrust_loading_horizontal_accel(wing_loading, deltamach, deltat, air_temp, dyn_press, CD0, K1, alpha, beta, CDR=0, K2=0, g0=9.81, gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Thrust loading during constant altitude and constant speed turn
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    """
    
    first  = K1*beta/dyn_press
    second = (CD0+CDR)/beta*dyn_press
    third = np.sqrt(gamma*R_gas*air_temp)*(deltamach)/g0/deltat
    thrust_loading = beta/alpha*(first*wing_loading + K2 + second/wing_loading + third) 

    return thrust_loading

def thrust_loading_landing(wing_loading, densitys, CLmax, CD, alpha, beta, friction_coeff, vel_ratio, distance, rotation_time=3, CL=0, CDR=0, g0=9.81):
    """
    Notes
    -----
    Thrust loading during landing free roll and braking
    
    Thrust loading = Thrust/Wto
    Wing Loading = Wto/S
    standard rotation time is 3 secs
    epsilon is the sum of the drags on the aircraft
    kto is the ratio of takeoff veloicty to stall velocity
    """
    epsilon = CD + CDR + friction_coeff*CL
    third = beta/densitys/g0/epsilon
    second = rotation_time*vel_ratio*np.sqrt(2*beta/densitys/CLmax)
    
    if alpha != 0:
        first = (distance-second*np.sqrt(wing_loading))/third/wing_loading
        thrust_loading = (epsilon*vel_ratio**2/(np.exp(first)-1)/CLmax-friction_coeff)*beta/alpha
    else:
        fourth = third*np.log(1+epsilon/friction_coeff/CLmax*vel_ratio**2)
        thrust_loading = ((-second+np.sqrt(second**2+4*fourth*distance))/2/fourth)**2

    return thrust_loading

def main():
    
    wing_loading1a = np.linspace(20, 120, 30)
    wing_loading1 = units.convert_pressure(wing_loading1a/144, 'Pa')
    
    densitys = 1.05498044 # kg/m^3
    # densitys =.002047
    CLmax = 2
    TR = np.array([1.07])
    beta1 = 1
    vel_ratio_to = 1.2
    takeoff_distance = 1500*12*.0254 # m
    mach1 = .1
    mach2 = 1.5
    mach3 = 1.6
    mach4 = .9
    CD = .3612
    mu = .05 # wet hard surface
    K12 = .27
    K13 = .288
    K14 = .18
    K15 = .216
    beta2 = .78
    load = 5
    dyn_press2 = units.convert_pressure(991.6/144, 'Pa')
    dyn_press3 = units.convert_pressure(1128/144, 'Pa')
    dyn_press4 = units.convert_pressure(357/144, 'Pa')
    dyn_press5 = units.convert_pressure(634.6/144, 'Pa')
    air_temp = get_atmos(30000)[0]*288.15
    CD02 = .028
    CD03 = .028
    CD04 = .016
    for i, tr in enumerate(TR):
        alpha1 = engine_thrust_lapse(tr, mach1, 1, 2, 1)
        alpha2 = engine_thrust_lapse(tr, mach2, 1, 30, 0)
        alpha3 = engine_thrust_lapse(tr, mach3, 1, 30, 1)
        alpha4 = engine_thrust_lapse(tr, mach4, 1, 30, 1)
        alpha5 = engine_thrust_lapse(tr, 1.2, 1, 30, 1)
        alpha7 = engine_thrust_lapse(tr, 1.8, 1, 40, 1)

        thrust_loading1 = thrust_loading_takeoff(wing_loading1, densitys, CLmax, CD, alpha1, beta1, mu, vel_ratio_to, takeoff_distance)
        thrust_loading2 = thrust_loading_cruise(wing_loading1, dyn_press2, CD02, K12, alpha2, beta2)
        thrust_loading3 = thrust_loading_turn(wing_loading1, load, dyn_press3, CD03, K13, alpha3, beta2)
        thrust_loading4 = thrust_loading_turn(wing_loading1, load, dyn_press4, CD04, K14, alpha4, beta2)
        thrust_loading5 = thrust_loading_horizontal_accel(wing_loading1, .8, 50, air_temp, dyn_press5, CD03, K15, alpha5, beta2)
        thrust_loading6 = thrust_loading_landing(wing_loading1a, .002047, 2, .8123, 0, .56, .18, 1.15, 1500, g0=32.17)
        thrust_loading7 = thrust_loading_cruise(wing_loading1a, 891.8, .028, .324, alpha7, .78)

        plt.plot(wing_loading1a, thrust_loading1, label="Takeoff", linestyle='--')
        plt.plot(wing_loading1a, thrust_loading2, label="Cruise", linestyle='--')
        plt.plot(wing_loading1a, thrust_loading3, label="Turn 1", linestyle='--')
        plt.plot(wing_loading1a, thrust_loading4, label="Turn 2", linestyle='--')
        plt.plot(wing_loading1a, thrust_loading5, label="Accel", linestyle='--')
        # plt.plot(wing_loading1a, thrust_loading6, label="Landing", linestyle='--')
        plt.vlines(thrust_loading6, .4, 1.6, label="Landing", linestyle='--')
        plt.plot(wing_loading1a, thrust_loading7, label="Max Mach", linestyle='--')

    plt.legend()
    plt.ylim([.4, 1.6])
    plt.xlim([20, 120])
    plt.show()

if __name__ == '__main__':
    main()