
import numpy as np
import USAtmos1976
import Unit_Conversions as units
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

Pstd = 101325 # Pa
Tstd = 288.15 # K
Dstd = 1.225 # kg/m^3
spdsnd_std = 340.3 # m/s

def getAtmosphere(altitude, day='std'):
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
            ValueError('Not a recognized day type. Use std, hot, trop, or cold')
    
    thetas = temp/288.15
    deltas = press/101325
    return (thetas, deltas, temp, press, density)

@dataclass
class GasTurbineEngine:
    thetat: int = None # total temperature : Tt/Tstd : theta_0
    thetas: int = None # static temperature : Ts/Tstd : theta
    deltat: int = None # total pressure : Pt/Pstd : delta_0
    deltas: int = None # static pressure : Ps/Pstd : delta
    sigmat: int = None # total density : densityt/Dstd : sigma_0
    sigmat: int = None # static density : Ds/Dstd : sigma
    throttle_ratio: int = None # throttle ratio Tt4max/Tt4SLS TR
    thrust_lapse: int = None # installed thrust lapse : thrust/sea_level_thrust
    gammac: int = None # average gamma for compression process
    __Rgas: int = 287.05
    theta_break: int = None # control system maximum Tt4 and compressor pressure ratio : theta_0_break
    Tt4max: int = 2400 # Rankine
    Tt4sls: int = 0

    engine_type: int = None
    engine_mode: int = None

    mach0: int = None
    theta0break: int = None
    """
    Variables 
    ---------
    variable name : actual name : quantity/description : in book

    wing_loading : wing loading : takeoff_weight/wing_planform_area : Wto/S
    thrust_load : thrust laoding : thrust/takeoff_weight : T/Wto
    engine_type : engine type : tj, lbtf, tp, etc : 
    engine_mode : engine mode : wet or dry : 
    
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
    velocity_ratio : velocity ratio : velocity/stall_velocity : k

    dynpress : dynamic pressure : density*velocity^2/2 : q
    
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

    def calcMFP(mach, gamma, Rgas):
        """
        Returns
        -------
        mfp : mass flow parameter m*sqrt(Tt)/Pt/area
        """
        exponent = (gamma+1)/2/(1-gamma)
        mfp = mach*np.sqrt(gamma/Rgas)*((1+(gamma-1)/2)*mach**2)**exponent
        return mfp

    def calcThetaBreak(self):
        """Theta break is the point at which the engine control system
        switches between limiting compression ratio and Tt4. It is the same
        think as throttle ratio"""
        
        theta_break = self.Tt4max/self.Tt4sls
        return theta_break

    def calcMach0Break(self):
        """Mach number at awhich theta0_break occurs"""
        self.mach0break = np.sqrt(2/(self.gammac-1)*(self.theta0break-1))
        return self.mach0break

    def calcThrustLapse(self, day, alt, gamma=1.4):
        """
        Thrust_lapse = thrust/thrust_sea_level
        Installed full throttle thrust lapse depends on altitude, speed, and
        whether or not an afterburner is operating. 

        page 38
        
        Parameters
        ----------
        engine_type : either of the following
            - hbtf - High bypass turbofan
            - lbtf - Low bypass turbofan
            - tj - Turbojet
            - tp - Turboprop
        engine_mode : either of the following
            - mil - military power
            - max - max power
        alt : altitude in feet
        """

        thetas = getAtmosphere(alt, day)[0]
        deltas = getAtmosphere(alt, day)[1]

        # normalized total temperature and total pressure
        norm_Tt = thetas*(1+(gamma-1)/2*self.mach0**2)
        norm_Pt = deltas*(1+(gamma-1)/2*self.mach0**2)**(gamma/(gamma-1))

        if norm_Tt <= self.throttle_ratio:
            zeta = 0
        else:
            zeta = 1

        match self.engine_type:
            case 'hbtf': # high bypass turbofan
                thrust_lapse = norm_Pt*(1-0.49*np.sqrt(self.mach0) - 3*(norm_Tt-self.throttle_ratio)/
                                        (1.5+self.mach0)*zeta)
            case 'lbtf': # low bypass turbofan
                if self.engine_mode == 'mil': # military power
                    thrust_lapse = 0.6*norm_Pt*(1-3.8*(norm_Tt-self.throttle_ratio)/norm_Tt*zeta)
                elif self.engine_mode == 'max': # maximum power
                    thrust_lapse = norm_Pt*(1-3.5*(norm_Tt-self.throttle_ratio)/norm_Tt*zeta)
            case 'tj': # turbojet
                if self.engine_mode == 'mil': # military power
                    thrust_lapse = 0.8*norm_Pt*(1-.16*np.sqrt(self.mach0) - 
                                                25/norm_Tt*(norm_Tt-self.throttle_ratio)
                                                /(9+self.mach0)*zeta)
                elif self.engine_mode == 'max': # maximum power
                    thrust_lapse = norm_Pt*(1-0.3*(norm_Tt-1)-0.1*np.sqrt(self.mach0) -
                                            1.5*(norm_Tt-self.throttle_ratio)/norm_Tt*zeta)
            case 'tp': # turboprop
                thrust_lapse = norm_Pt # less than mach = 0.1
                thrust_lapse = norm_Pt*(1-0.96*(self.mach0-1)**0.25 - 3/8.13*(
                    norm_Tt-self.throttle_ratio)/(self.mach0-0.1)*zeta)
        
        return thrust_lapse

    def getInitialTSFC(self, TsTstd):
        """
        Provide an initial estimate for the thrust specific fuel
        consumption units of 1/hour

        Parameters
        ----------
        TsTstd : Ts/Tstd

        page 71
        """
        match self.engine_type:
            case 'hbtf': # high bypass turbofan
                tsfc = (0.45+0.54*self.mach0)*np.sqrt(TsTstd)
            case 'lbtf': # low bypass turbofan
                if self.engine_mode == 'mil': # military power
                    tsfc = (0.9+0.30*self.mach0)*np.sqrt(TsTstd)
                elif self.engine_mode == 'max': # maximum power
                    tsfc = (1.6+0.27*self.mach0)*np.sqrt(TsTstd)
            case 'tj': # turbojet
                if self.engine_mode == 'mil': # military power
                    tsfc = (1.1+0.30*self.mach0)*np.sqrt(TsTstd)
                elif self.engine_mode == 'max': # maximum power
                    tsfc = (1.5+0.23*self.mach0)*np.sqrt(TsTstd)
            case 'tp': # turboprop
                tsfc = (0.18+0.8*self.mach0)*np.sqrt(TsTstd)

        return tsfc

    def getEngineConstants(self):
        """
        Get engine specific constants for mission analysis

        page 71

        Parameters
        ----------
        engine_type : either of the following
            - hbtf - high bypass turbofan constants returned
            - lbtf - high bypass turbofan constants returned
            - tj - turbojet constants returned
            - tp - turboprop constants returned
        engine_mode : either of the following
            - mil - military power constants returned
            - max - maximum power constants returned

        Returns
        -------
        0: C1 units 1/hour
        1: C2 units 1/hour
        """
        match self.engine_type:
            case 'hbtf': # high bypass turbofan
                consts = [0.45, 0.54]
            case 'lbtf': # low bypass turbofan
                if self.engine_mode == 'mil': # military power
                    consts = [0.9, 0.30]
                elif self.engine_mode == 'max': # maximum power
                    consts = [1.6, 0.27]
            case 'tj': # turbojet
                if self.engine_mode == 'mil': # military power
                    consts = [1.1, 0.30]
                elif self.engine_mode == 'max': # maximum power
                    consts = [1.5, 0.23]
            case 'tp': # turboprop
                consts = [0.18, 0.8]
        
        return consts
    
@dataclass
class AerospaceVehicle:
    # constants
    gamma: int = 1.4
    mach0: int = None
    altitude: int = None
    
    # aerodynamic attributes
    CDR: int = 0
    CD0: int = None
    K1: int = None
    CD: int = None # coefficient of drag
    CLmax: int = None
    engine: GasTurbineEngine = None

    # chaning attributes
    instwf: int = None # instantaneous weight fraction, weight/takeoff weight
    thrust_load: int = None # sea level thrust to takeoff weight ratio 
    wing_loading: int = None # takeoff weight to planform area ratio
    takeoff_distance: int = None

    # other attributes
    mach_crit: int = 0.9 # mach number at which divergence drag starts
    
    def calcDynamicPressure(self, alt, mach, gamma=1.4):
        """
        Notes
        -----
        Altitude in feet. Can be a nummpy array
        
        page 11
        
        Returns
        -------
        dynamic pressure in Pascals

        Parameters
        ----------
        alt : altitude in feet
        """

        deltas = self.getAtmosphere(alt, 'std')[1]
        return deltas*101325*mach**2*gamma/2
    
    def calcDragPolar(self, velocity_ratio=1.25, K2=0, CL=0, suppress=False):
        """
        Notes
        -----
        Calculate drag polar
        K1, K2, CD0 are assumed known
        high performance aircraft CLmin~=0 K2~=0
        K2 = 0 for uncambered
        0.001 <= K" <= 0.03
        0.1 <= CLmin <= 0.3
        7 <= AR <= 10
        0.75 <= efficiency <= 0.85

        Lift-drag polar breakdown
        CD = CDmin + K'CL^2 + K"(CL - CLmin)^2
        CD = K1CL^2 + K2CL + CDO
        
        page 25,37

        Parameters
        ----------
        self.mach0 : freestream mach number
        CLmax : maximum coefficient of lift
        velocity_ratio : V = velocity_ratio*Vstall

        Returns
        ------
        0: CD0
        1: K1
        2: CD    
        """

        if 0 <= self.mach0 <= 1:
            self.K1 = 0.18
        elif 1 < self.mach0 <= 2:
            self.K1 = (.18 - .36)/(1-2)*(self.mach0-1)+.18
        else:
            self.K1 = 0
            if not suppress:
                ValueError("Flight mach number not in range for turbine engines.")

        if 0 <= self.mach0 <= 0.8:
            self.CD0 = .014
        elif  0.8 < self.mach0 <= 0.9:
            self.CD0 = (.014 - .016)/(0.8-.9)*(self.mach0-0.8)+.014
        elif 0.9 < self.mach0 <= 1.1:
            self.CD0 = (.016 - .026)/(0.9-1.1)*(self.mach0-0.9)+.016
        elif 1.1 < self.mach0 <= 1.2:
            self.CD0 = (.026 - .028)/(1.1-1.2)*(self.mach0-1.1)+.026
        elif 1.2 < self.mach0 <= 2:
            self.CD0 = .028
        else:
            self.CD0 = 0
            if not suppress:
                ValueError("Flight mach number not in range for turbine engines.")
            
        if not CL: CL = self.CLmax/velocity_ratio**2
        self.CD = self.K1*CL**2+K2*CL+self.CD0

        return (self.CD0, self.K1, self.CD)
    
    def calcEmptyWeightFraction(takeoff_weight, type='fighter'):
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
                empty_wghtfrac = 2.34*takeoff_weight**-0.13
            case 'cargo':
                empty_wghtfrac = 1.26*takeoff_weight**-0.08
            case 'passenger':
                empty_wghtfrac = 1.02*takeoff_weight**-0.06
            case 'twin prop':
                empty_wghtfrac = 0.96*takeoff_weight**-0.05

        return empty_wghtfrac
    
    def calcTakeoffVelocity(self, Ds, kto=1.2):
        """
        Notes
        -----
        Calculates the takeoff velocity
        
        page 29
        """
        vel_to = kto*np.sqrt(2*self.instwf/Ds/self.CLmax*self.wing_loading)
        print(284, kto, self.instwf, Ds, self.CLmax, self.wing_loading)
        return vel_to

    def calcThrustLoadingTurn(self, load_factor, alt, mach, CDR=0, K2=0):
        """
        Notes
        -----
        Thrust loading during constant altitude and constant speed turn
        
        Thrust loading = Thrust/Wto
        Wing Loading = Wto/S

        page 26
        """
        dyn_press = self.calcDynamicPressure(alt, mach)
        CD0, K1, CD = self.calcDragPolar(mach, self.CLmax)

        help1  = K1*load_factor**2*self.instwf/dyn_press
        help2 = (CD0+CDR)/self.instwf*dyn_press
        thrust_load = self.instwf/self.engine.thrust_lapse*(help1*self.wing_loading + K2*load_factor + 
                                    help2/self.wing_loading) 

        return thrust_load

    def calcExcessPower(self, load_factor, velocity, alt, mach,  K2=0, CDR=0):
        """
        Notes
        -----
        Calculates excess power for the aircraft at the defined thrust
        loading, wing loading, throttle ratio for any instantenous weight fraction, load factor
        across the flight envelope (altitude vs altitude)

        page 48

        Returns
        -------
        excess power

        Parameters
        ----------
        alt : altitude in feet no regard for other units
        all other parameters units in SI
        """
        dyn_press = self.calcDynamicPressure(alt, mach)
        CD0, K1, CD = self.calcDragPolar(mach, CLmax, suppress=True)

        help1  = self.engine.thrust_lapse/self.instwf*self.thrust_load
        help2 = -K1*load_factor**2*instwf/dyn_press*wing_loading
        help3  = -K2*load_factor
        fourth = -(CD0+CDR)/instwf/wing_loading*dyn_press
        excess_power = velocity*(help1+help2+help3+fourth)

        return excess_power

@dataclass
class MissionAnalysis:
    """Mission description on page 13
    Thrust Loading equations are on page 25"""
    mission_phases: list = None
    vehicle: AerospaceVehicle = None
    kto: int = None # velocity ratio at takeoff
    ktd: int = None # velocity ratio at touchdown
    g0: float = 9.81 # m/s^2

    def performContstraintAnalysis(self, engine_type):
        """
        Notes
        -----
        Need to add remining constraint equations form chapter 2 but
        currently not necessary for remainder of book walkthrough
        """

        # Mission Phase 1: Takeoff, no obstacle
        
        wing_loading1a = np.linspace(20, 120, 30)
        wing_loading1 = units.convert_pressure(wing_loading1a/144, 'Pa')
        
        Ds = 1.05498044 # kg/m^3
        TR = 1.07
        instwf = [1, .78, .78, .78, .78, .56, .78]
        takeoff_distance = 1500*12*.0254 # m
        mach = [.1, 1.5, 1.6, .9, 1.2, 0, 1.8]
        alt = [2000, 30000, 30000, 30000, 30000, 0, 40000]

        load_turn1 = load_turn2 = 5 # g

        deltamach = .8
        deltat = 50

        air_temp = self.getAtmosphere(30000, 'std')[0]*288.15
        kto = 1.2
        ktd = 1.15
        muto = .05
        mutd = .18
        CLmax = 2

        landing_drag = .8123
        thrust_reverse = 0

        thrust_lapse1 = self.engine.calcThrustLapse(TR, mach[0], engine_type, 'max', alt[0])
        thrust_lapse2 = self.engine.calcThrustLapse(TR, mach[1], engine_type, 'mil', alt[1])
        thrust_lapse3 = self.engine.calcThrustLapse(TR, mach[2], engine_type, 'max', alt[2])
        thrust_lapse4 = self.engine.calcThrustLapse(TR, mach[3], engine_type, 'max', alt[3])
        thrust_lapse5 = self.engine.calcThrustLapse(TR, mach[4], engine_type, 'max', alt[4])
        thrust_lapse7 = self.engine.calcThrustLapse(TR, mach[6], engine_type, 'max', alt[6])
        # All thrust loading equations have been validated
        thrust_load1 = self.calcThrustLoadingTO(wing_loading1, Ds, CLmax, mach[0], thrust_lapse1, 
                                    instwf[0], muto, kto, takeoff_distance) 
        thrust_load2 = self.calcThrustLoadingCruise(wing_loading1, alt[1], mach[1], CLmax, thrust_lapse2,
                                    instwf[1]) 
        thrust_load3 = self.calcThrustLoadingTurn(wing_loading1, load_turn1, alt[2], mach[2], 
                                CLmax, thrust_lapse3, instwf[2]) 
        thrust_load4 = self.calcThrustLoadingTurn(wing_loading1, load_turn2, alt[3], mach[3], CLmax, 
                                thrust_lapse4, instwf[3]) 
        thrust_load5 = self.calcThrustLoadingAccel(wing_loading1, deltamach, deltat, air_temp,
                                        alt[4], mach[4], CLmax, thrust_lapse5, instwf[4]) 
        thrust_load6 = self.calcThrustLoadingLanding(wing_loading1, Ds, CLmax, landing_drag, 
                                    thrust_reverse, instwf[5], mutd, ktd, 
                                    takeoff_distance)
        thrust_load7 = self.calcThrustLoadingCruise(wing_loading1, alt[6], mach[6], CLmax, thrust_lapse7, 
                                    instwf[6]) 
        
        thrust_load6 = units.convert_pressure(thrust_load6, "psi")*144

        plt.plot(wing_loading1a, thrust_load1, label="Takeoff", linestyle='--')
        plt.plot(wing_loading1a, thrust_load2, label="Cruise", linestyle='--')
        plt.plot(wing_loading1a, thrust_load3, label="Turn 1", linestyle='--')
        plt.plot(wing_loading1a, thrust_load4, label="Turn 2", linestyle='--')
        plt.plot(wing_loading1a, thrust_load5, label="Accel", linestyle='--')
        plt.vlines(thrust_load6, .4, 1.6, label="Landing", linestyle='--')
        plt.plot(wing_loading1a, thrust_load7, label="Max Mach", linestyle='--')

        plt.legend()
        plt.ylim([.4, 1.6])
        plt.xlim([20, 120])
        plt.show()
        plt.close()

    def performPowerAnalysis(self, throttle_ratio, engine_type, CLmax):
        """
        Notes
        -----
        Function to compare this codes results to book results
        
        Mimics plot on page 48 but slight indent on righgt corner of envelope
        Single point calculation performed after plotting gives 316 ft/s
        book gives 320 ft/s

        More resolution needed but takes to long to render plot, convert to
        C or Fortran
        """

        load_factor = 1
        instwf = .97
        wing_loading = units.convert_pressure(wing_loading/144, 'Pa')
        N = 50
        velocity = np.linspace(100, 1900, N)*.0254*12 # ft/sec
        altitude = np.linspace(0, 60000, N) # ft
        power = np.zeros((N, N))
        
        for i, alt in enumerate(altitude):
            for j, vel in enumerate(velocity):
                mach = (vel/USAtmos1976.stdday(alt).speed_of_sound)[0]
                thrust_lapse = self.calcThrustLapse(self.throttle_ratio, mach, engine_type, 'mil', alt)
                power[i][j] = self.calcExcessPower(wing_loading, thrust_load, load_factor, vel, alt, 
                                        mach, CLmax, thrust_lapse, instwf)
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
        thrust_lapse = self.calcThrustLapse(self.throttle_ratio, mach, engine_type, 'mil', alt)
        pwr = calcExcessPower(wing_loading, thrust_load, load_factor, vel, alt, mach, 
                        CLmax, thrust_lapse, instwf)
        print(pwr/.0254/12)

    def performMissionAnalysis(self, thrust_load, wing_loading, throttle_ratio, engine_type, CLmax):
        """
        Notes
        -----
        Perform a weight analysis over the entire mission set
        """
        wing_loading = units.convert_pressure(wing_loading/144, 'Pa')
        alt_to = 2000
        instwfs_to = self.getWeightFractionTakeoff(alt_to, thrust_load, wing_loading, self.throttle_ratio, engine_type, CLmax)
        instwf_to = instwfs_to[0]*instwfs_to[1]*instwfs_to[2]
        print(instwf_to)
        mach = .4410
        veli = 211.1*.0254*12
        velf = 812*.0254*12
        alt = 2000
        thrust_lapse = self.calcThrustLapse(self.throttle_ratio, mach, engine_type, 'mil', 'trop', 2000)
        self.getWeightFractionAccel(instwf_to, alt, wing_loading, mach, thrust_lapse, CLmax, velf, veli, engine_type, 'mil')

@dataclass
class WarmupTakeoff(MissionAnalysis):
    """
    Warm-up and takeoff, field is at 2000 ft pressure altitude (PA) with air
    temperature of 100°E Fuel allowance is 5 min at idle power for taxi and 1 rain
    at military power (mil power) for warm-up. Takeoff ground roll plus 3 s
    rotation distance must be < 1500 ft on wet, hard surface runway (mu = 0.05),
    Vro = 1.2 Vstall"""
    friction_coeff: float = 0.05
    altitude: float = 2000 # ft
    mach: float = 0.1
    air_temp = 100 # F
    takeoff_distance: float = 1500 # ft
    velocity_ratio: float = 1.2 # Vro/Vstall
    weight_fraction: float = None # weight fraction at takeoff
    inst_weightfrac: float = None # instantaneous weight fraction
    rotation_time: float = 3 # seconds # time to for vehicle to rotate when taking off
    wing_loading:float = None

    def calcThrustLoading(self, Ds):
        """
        Thrust loading during takeoff acceleration up to rotation
        instantaneous weight is given by W = instwf*WTO.
        Instantaneous weight fraction, depends on how much fuel has been
        consumed and payload delivered; more AED2e p24

        Thrust loading = Thrust/Wto
        Wing Loading = Wto/S
        standard rotation time is 3 secs
        epsilon is the sum of the drags on the aircraft
        kto is the ratio of takeoff veloicty to stall velocity

        page 34
        """
        CD0, K1, CD = self.vehicle.calcDragPolar(self.mach, self.vehicle.CLmax, self.velocity_ratio)
        CL = self.vehicle.CLmax/self.velocity_ratio**2
        CDR = self.friction_coeff*CL # coefficient of additional drags
        takeoff_resistance = CD + CDR + self.friction_coeff*CL # xi_to

        a0 = self.inst_weightfrac/Ds/self.g0/takeoff_resistance
        b = self.rotation_time*self.velocity_ratio*np.sqrt(2*self.inst_weightfrac/Ds/self.vehicle.CLmax)
        help1 = -(self.takeoff_distance-b*np.sqrt(self.wing_loading))/a0/self.wing_loading
        
        thrust_load = (takeoff_resistance*self.velocity_ratio**2/(1-np.exp(help1))/self.vehicle.CLmax+
                        self.friction_coeff)*self.inst_weightfrac/self.vehicle.engine.thrust_lapse

        return thrust_load

    def calcWeightFractionWarmup(self, deltat, thrust_load, thrust_lapse):
        """
        Notes
        -----
        Calculates weight fraction of aircraft after warm-up

        page 67
        """
        thetas = self.getAtmosphere(self.altitude)[0]
        C1, _ = self.vehicle.engine.getEngineConstants()

        wght_frac = 1 - C1/3600*np.sqrt(thetas)*thrust_lapse/self.instwf*thrust_load*deltat
        
        return wght_frac
    
    def calcWeightFractionAccel(self, drgcf, frccf, mach, mach_to, vel,
                             g0=9.81):
        """
        Notes
        -----
        Calculates weight fraction of aircraft after takeoff acceleration

        page 63
        """
        thetas = self.getAtmosphere(self.altitude)[0]
        C1, C2 = self.aircraft.engine.getEngineConstants()
        dyn_press = self.calcDynamicPressure(self.altitude, mach_to)

        ttldrag_thrst = (drgcf*dyn_press/self.instwf/self.wing_loading+frccf)*(self.instwf/self.engine.thrust_lapse/self.thrust_load)
        print(700, self.instwf, self.engine.thrust_lapse, self.thrust_load)
        weight_frac = np.exp(-(C1+C2*mach)/3600*np.sqrt(thetas)/g0*vel/
                            (1-ttldrag_thrst))
        return weight_frac
    
    def getWeightFractionRotation(self, mach_to, rot_time=3):
        """
        Notes
        -----
        Calculates weight fraction of aircraft after takeoff rotation

        page 68
        """
        thetas = self.getAtmosphere(self.altitude)[0]
        C1, C2 = self.aircraft.engine.getEngineConstants()

        wght_frac = 1-(C1+C2*mach_to)/3600*np.sqrt(thetas)*self.engine.thrust_lapse/self.instwf*(self.thrust_load*rot_time)

        return wght_frac

    def getWeightFractionAccel(self, instwf_prev, alt, mach, thrust_lapse, velf, veli, CDR=0, g0=9.81, gamma=1.4):
        """
        Notes
        -----
        Calculates weight fraction of aircraft after horiztonal acceleration

        p 62
        """

        CL = 2*instwf_prev*(self.wing_loading/gamma/Pstd/self.getAtmosphere(alt, 'std')[1]/mach**2)
        CD0, K1, CD = self.calcDragPolar(mach, self.CLmax, CL=CL)
        C1, C2 = self.engine.getEngineConstants()
        # print(CL, C1/mach+C2, CD)
        print(thrust_lapse, instwf_prev)
        help1 = -(C1/mach+C2)*3600/spdsnd_std
        dvel = velf**2 - veli**2
        help3 = (CD+CDR)/CL/(instwf_prev/thrust_lapse)/self.wing_loading
        # print(help1/3600*spdsnd_std)
        # print(dvel/2/g0/12/.0254)
        # print(help3, CD, CL, instwf_prev, thrust_lapse)
        help2 = dvel/2/g0/(1-help3)
        instwf = np.exp(help1*help2)
        print(instwf)
        return instwf
    
    def getWeightFraction(self, alt_to):
        """
        Notes
        -----
        Calculates the weight fraction the takeoff sequence: warm-up, acceleration, and rotation
        Returns three weight fraction values: warm-up, acceleration, and rotation
        """
        
        mach_wu = 0
        instwf_wu = 1
        deltat = 60 # s
        engine_mode_wu = 'mil'
        thrust_lapse_wu = self.engine.calcThrustLapse('std', alt_to)
        instwf_wu = self.calcWeightFractionWarmup(deltat, engine_mode_wu, thrust_lapse_wu, instwf_wu)

        mach_to = .1819
        mach_toa = .1
        drgcf_toa = .3612
        frccf_toa = .05
        engine_mode_toa = 'max'
        Ds = self.getAtmosphere(alt_to, 'std')[4]
        vel_to = self.calcTakeoffVelocity(instwf_wu, Ds)
        thrust_lapse_toa = self.calcThrustLapse(self.engine.throttle_ratio, mach_toa, engine_mode_toa, 'std', alt_to)
        instwf_toa = self.calcWeightFractionTOAccel(alt_to, drgcf_toa, frccf_toa, mach_toa, mach_to, vel_to,
                        engine_mode_toa, thrust_lapse_toa, instwf_wu)
        
        instwf1 = instwf_toa*instwf_wu
        thrust_lapse_torot = self.calcThrustLapse(self.engine.throttle_ratio, mach_to, engine_mode_toa, 'std', alt_to)
        instwf_torot = self.getWeightFractionTORotation(alt_to, mach_to, engine_mode_toa, thrust_lapse_torot, instwf1)
        
        return (instwf_wu, instwf_toa, instwf_torot)

@dataclass
class AccelerateClimb(MissionAnalysis):
    """"
    Accelerate to climb speed and perform a minimum time climb in mil power to
    best cruise Mach number and best cruise altitude conditions (BCM/BCA)"""

    starting_mach: int = None
    ending_mach: int = None
    altitude: int = None # ft
    mach: int = None
    
    def calcThrustLoadingAccel(self, deltamach, deltat, air_temp, instwf, 
                            CDR=0, K2=0, g0=9.81, gamma=1.4, Rgas=287.05):
        """
        Thrust loading during acceleration before climb
        
        page 27
        """
        dyn_press = self.aircraft.calcDynamicPressure(self.altitude, self.mach)
        CD0, K1, CD = self.aircraft.calcDragPolar(self.mach, self.aircraft.CLmax)

        help1  = K1*instwf/dyn_press
        help2 = (CD0+CDR)/instwf*dyn_press
        help3 = np.sqrt(gamma*Rgas*air_temp)*(deltamach)/g0/deltat
        thrust_load = instwf/thrust_lapse*(help1*wing_loading + K2 + 
                                    help2/wing_loading + help3) 

        return thrust_load

@dataclass
class SupersonicCruiseClimb(MissionAnalysis):
    """
    Subsonic cruise climb at BCM/BCA until total range for climb and cruise climb
    is 150 n miles"""
    """
    Subsonic cruise climb at BCM/BCA until total range from the end of combat
    equals 150 n miles. 
    """
    range: float = 150 # n miles

    def getBCM_BCA(self):
        """
        Notes
        -----
        Obtain the best cruising altitude for subsonic speeds. The best
        cruise mach number is M=0.9 after which divergence drag begins,
        for fighter aircraft
        """
        CD0, K1, _ = self.calcDragPolar(self.mach_crit, self.CLmax)
        bca = 2*self.instwf*self.wing_loading/self.gamma/self.Pstd/self.mach_crit*2*np.sqrt((CD0+self.CDR)/K1)

        return bca

@dataclass
class Descent(MissionAnalysis):
    """
    Descend to 30,000 ft. No range/fuel/time credit for descent."""
    """
    Descend to 10,000 ft. No time/fuel/distance credit"""
    pass

@dataclass
class CombatAirPatrol(MissionAnalysis):
    """
    Perform a combat air patrol (CAP) loiter for 20 min at 30,000 ft and Mach
    number for best endurance."""
    pass

@dataclass
class SupersonicPenetration(MissionAnalysis):
    """
    Supersonic penetration at 30,000 ft and M = 1.5 to combat arena. Range =
    100 n miles. Penetration should be done as a military power (i.e., no
    afterburning) supercruise if possible. """
    altitude: float = 30000 # ft
    mach: float = 1.5 # M
    engine_mode: str = 'mil' # military power
    range: float = 100 # n miles

    def calcThrustLoadingCruise(self, alt, mach, CDR=0, K2=0):
        """
        Thrust loading during cruise

        page 26
        """
        dyn_press = self.calcDynamicPressure(alt, mach)
        CD0, K1, CD = self.calcDragPolar(mach, self.CLmax)
        thrust_load = self.instwf/self.engine.thrust_lapse*(K1*self.instwf/dyn_press*self.wing_loading + K2 + (CD0+CDR)/self.instwf*
                            dyn_press/self.wing_loading) 

        return thrust_load

@dataclass
class Combat(MissionAnalysis):
    """
    Combat is modeled by the following:
    Fire 2 AMRAAMs
    Perform one 360 deg, 5g sustained turn at 30,000 ft. M = 1.60
    Perform two 360 deg, 5g sustained turns at 30,000 ft. M = 0.90
    Accelerate from M = 0.80 to M = 1.60 at 30,000 ft in maximum power
    (max power)
    Fire 2 AIM-9Ls and 1/2 of ammunition
    No range credit is given for combat maneuvers.
    Conditions at end of combat are M = 1.5 at 30,000 ft."""
    altitude: float = 30000 # ft
    mach: float = 1.6 # M
    engine_mode: str = 'max' # maximum power

@dataclass
class EscapeDash(MissionAnalysis):
    """
    Escape dash, at M = 1.5 and 30,000 ft for 25 n miles. Dash should be done as a
    mil power supercruise if possible. 
    """
    altitude: float = 30000 # ft
    mach: float = 1.5 # M
    range: float = 25 # n miles
    engine_mode: str = 'mil' # military power

    def calcThrustLoadingAccel(wing_loading, deltamach, deltat, air_temp, alt, mach, 
                            CLmax, thrust_lapse, instwf, 
                            CDR=0, K2=0, g0=9.81, gamma=1.4, Rgas=287.05):
        """
        Notes
        -----
        Thrust loading during accletation
        
        Thrust loading = Thrust/Wto
        Wing Loading = Wto/S

        page 27
        """
        dyn_press = calcDynamicPressure(alt, mach)
        CD0, K1, CD = calcDragPolar(mach, CLmax)

        help1  = K1*instwf/dyn_press
        help2 = (CD0+CDR)/instwf*dyn_press
        help3 = np.sqrt(gamma*Rgas*air_temp)*(deltamach)/g0/deltat
        thrust_load = instwf/thrust_lapse*(help1*wing_loading + K2 + 
                                    help2/wing_loading + help3) 

        return thrust_load

@dataclass
class Climb(MissionAnalysis):
    """
    Using mil power, perform a minimum time climb to BCM/BCA. (If the initial
    energy height exceeds the final, a constant energy height maneuver may be
    used. No distance credit for the climb.)"""
    pass

@dataclass
class Loiter(MissionAnalysis):
    """
    Loiter 20 min at 10,000 ft and Mach number for best endurance"""
    pass

@dataclass
class Landing(MissionAnalysis):
    """
    Descend and land, field is at 2000 ft PA, air temperature is 100°F. A 3 s free
    roll plus braking distance must be < 1500 ft. On wet, hard surface runway
    (/z8 = 0.18), VTD =1.15 VSTALL. """
    landing_distance: float = 1500 # ft
    altitude: float = 2000 # ft
    air_temp: float = 100 # F
    velocity_ratio: float = 1.15 # VTD/VSstall
    friction_coeff: float = 0.18 # friction coefficient

    def calcThrustLoadingLanding(wing_loading, Ds, CLmax, CD, thrust_lapse, 
                                instwf, frccf_ln, velocity_ratio, distance, 
                                rotation_time=3, CL=0, CDR=0, g0=9.81):
        """
        Notes
        -----
        Thrust loading during landing free roll and braking
        
        Thrust loading = Thrust/Wto
        Wing Loading = Wto/S
        standard rotation time is 3 secs
        epsilon is the sum of the drags on the aircraft
        kto is the ratio of takeoff veloicty to stall velocity

        page 31
        """
        epsilon = CD + CDR + frccf_ln*CL
        help3 = instwf/Ds/g0/epsilon
        help2 = rotation_time*velocity_ratio*np.sqrt(2*instwf/Ds/CLmax)
        
        if thrust_lapse != 0:
            help1 = (distance-help2*np.sqrt(wing_loading))/help3/wing_loading
            thrust_load = (epsilon*velocity_ratio**2/(np.exp(help1)-1)/CLmax-
                            frccf_ln)*instwf/thrust_lapse
        else:
            fourth = help3*np.log(1+epsilon/frccf_ln/CLmax*velocity_ratio**2)
            thrust_load = ((-help2+np.sqrt(help2**2+4*fourth*distance))/
                            2/fourth)**2

        return thrust_load
        
def main():
    AAF = AerospaceVehicle()
    engine_type = 'lbtf' # low bypass turbofan

    # performContstraintAnalysis(engine_type)

    # Aircraft intrinsic properties
    CLmax = 2
    # Chosen design point
    thrust_load = 1.25
    wing_loading = 64 # lbf/ft^2
    self.throttle_ratio = 1.07

    # power_anlys(thrust_load, wing_loading, self.throttle_ratio, engine_type, CLmax)
    performMissionAnalysis(thrust_load, wing_loading, self.throttle_ratio, engine_type, CLmax)

if __name__ == '__main__':
    # main()

    # thrust_lapse = calcThrustLapse(1.07, .4410, 'lbtf', 'mil', 'std', 2000)
    # print(thrust_lapse)
    hs = np.linspace(0, 30000, 4)/.0254/12
    t1 = getAtmosphere(hs, 'std')
    # t2 = getAtmosphere(30000, 'cold')
    # t3 = getAtmosphere(30000, 'hot')
    # t4 = getAtmosphere(30000, 'trop')

    # print(t1, t2, t3, t4)
    print(hs*.0254*12)
    print(t1)