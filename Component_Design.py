
import numpy as np
import IsentropicFlow as isenf
import RootFinding as rootfind
import Unit_Conversions as units

def _find_mach(mach, massflow, densityt, Tt, area, gamma, R_gas):
    """
    Notes
    -----
    Find Mach number from given massflow, total temperature, density, 
    and flow area. Mach number is the variable that is changing to 
    converge get the areas to match. Function to be used with the 
    newton2 solver in the RootFinding module. Function iterates even
    though only one unknown is present because rearraning the equations
    is too complicated
    
    Returns
    -------
    0: difference between flow area calculated from convservation of 
    mass and the given flow area
    
    Parameters
    ----------
    mach : mach number
    massflow : mass flow rate
    densityt : total density
    Tt : total temperature
    area : flow area to be matched

    Assumptions
    -----------
    0: 
    """
    densitys = isenf.static_density(mach, densityt, gamma)
    Ts = isenf.static_temperature(mach, Tt, gamma)
    sound_speed = np.sqrt(gamma*R_gas*Ts)
    return (massflow/densitys/mach/sound_speed - area,)

def _find_mach2(mach, massflow, densityt, Tt, flow_area, velocity_comp, gamma,
                 R_gas):
    """
    Notes
    -----
    Find Mach number from given massflow, total temperature, density, 
    and flow area and a given velocity component. Mach number is the 
    variable that is changing to converge get the velocitys to match. 
    The function is used when a certain velocity needs to be known to 
    satisfy mass conservation but there is a velocity component that 
    does not contribute to conservation of mass. As in when only the 
    axial component of exhaust velocity contributes to mass flow, and 
    not the azimthual component, in a turbine. The squared velocities 
    are compared because during iterations velocity component might, 
    impossibly, be greater than the total velocity of the exhaust.
    Function to be used with the newton2 solver in the RootFinding 
    module. Function iterates even though only one unknown is present 
    because rearraning the equations is too complicated
    
    Returns
    -------
    0: difference between squared velocities calculated from mach number
    equations and from mass flow equations
    1: velocity component
    
    Parameters
    ----------
    mach : mach number
    massflow : mass flow rate
    densityt : total density
    Tt : total temperature
    area : flow area to be matched
    velocity_component : velocity component that does not contirubte to
    mass flow

    Assumptions
    -----------
    0: 
    """
    Ts1 = isenf.static_temperature(mach, Tt, gamma)
    densitys1 = isenf.static_density(mach, densityt, gamma)
    cvelz1sq = (mach*np.sqrt(gamma*R_gas*Ts1))**2 - velocity_comp**2
    veclocity1zsq2 = (massflow/densitys1/flow_area)**2
    return (veclocity1zsq2 - cvelz1sq, np.sqrt(cvelz1sq))

def crit_pressure_coeff(mach0, gamma=1.4):
    """
    Notes
    -----
    Returns the pressure coefficient that would occur when the incoming
    flow is accelerated to sonic conditions
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
    To be used with newton2 rootfind to determine the mach number to
    which the flow would be accelerated if it experienced an area change
    from A0 to A1. Can be reworked to be solved faster and remove
    dependence on rootfind

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

def additive_drag(mach0:float, mach1:float, gamma=1.4):
    """
    Notes
    -----
    Calculates the inlet additive drag due to air intake obstructing the
    flow

    Returns
    -------
    drag

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

def prediff_mach(mach2, mach1, Pt31, Ptin, area31, dome_area, gamma):
    """
    Notes
    -----
    To be used with rootfind, newton2
    p517
    """
    first = 1 + (gamma-1)/2*mach1**2
    second = 1 + (gamma-1)/2*mach2**2
    mach22 = Pt31/Ptin*area31/dome_area*mach1*(second/first)**(
        (gamma+1)/(gamma-1)/2)

    return (mach22-mach2, mach22)

def _combustor_casing_calcs(Tt31, Pt31, airflow, ref_vel, pitch_diam, R_gas):
    """Calculations for determining casing values used in combustor
    function"""

    rhot31 = Pt31/R_gas/Tt31
    ref_area = airflow/rhot31/ref_vel
    ref_height = ref_area/np.pi/pitch_diam
    diam_inner_casing = (pitch_diam/2 - ref_height/2)*2
    diam_outer_casing = (pitch_diam/2 + ref_height/2)*2

    return (diam_inner_casing, diam_outer_casing, rhot31, ref_height, ref_area)

def _combustor_dome_calcs(flow_split, airflow, passage_vel, diam_inner_casing, 
                          diam_outer_casing, rhot31, ref_height, pitch_diam, 
                          dome_height_set=None):
    """Calculations for determining combustor dome values used in
    combustor function"""

    area_passage = flow_split*airflow/rhot31/passage_vel  
    diam_inner_pass = 2*np.sqrt(diam_inner_casing**2/4 + area_passage/np.pi)
    diam_outer_pass = 2*np.sqrt(diam_outer_casing**2/4 - area_passage/np.pi)
    height_inner_pass = (diam_inner_pass - diam_inner_casing)/2
    height_outer_pass = (diam_outer_casing - diam_outer_pass)/2

    dome_height = ref_height - height_inner_pass - height_outer_pass
    if dome_height_set: dome_height = dome_height_set
    dome_area = np.pi*dome_height*pitch_diam
    dome_vel = flow_split*airflow/rhot31/dome_area
    
    return (dome_height, dome_vel, diam_inner_pass, diam_outer_pass)

def _space_rate_calcs(Pt31, airflow, fuelflow, pitch_diam, LHV, 
                      comblen_domeheight, height_turbine_inlet, dome_height):
    """Calculations for determining space rate used in combustor
    function"""
        
    comb_length = comblen_domeheight*dome_height
    area_entr = dome_height
    area_exit = (dome_height+height_turbine_inlet)/2
    comb_vol = (area_entr+area_exit)*np.pi*pitch_diam*comb_length/2
    # print(comblen_domeheight, comb_length/.0254, comb_vol/(.0254*12)**3)

    fuel_air = fuelflow/airflow
    airflow_lb = units.convert_mass(airflow, 'lbm')
    Pt31_atm = Pt31/101300
    comb_vol_ft3 = comb_vol/(.0254*12)**3
    space_rate = 3600*fuel_air*airflow_lb*LHV/Pt31_atm/comb_vol_ft3
    # print('space rate', space_rate)

    return (space_rate, comb_length)
            
def _find_cz(velc1, Ts0, Pt1, w0, velct1, massflow, flowarea1, spool_speed1, 
             spool_speed2, cp_gas, gminusg, R_gas):
    """
    Notes
    -----
    The axial velocity of the air after the rotor is needed and must
    satisfy mass flow rates while also having a component in the
    azimuthal direction velc1 is the variable that is changing to
    converge the velocities Follows the same numbering as parent
    function compressor_blade_design

    Returns
    -------
    0: difference of axial velcoity guessed and calculated from mass
    flow equations
    1: total temperature at station 1
    2: static temperature at station 1

    Parameters
    ----------
    velc1 : total absolute velocity at after rotor
    Ts0 : static temperature before rotor
    Pt1 : total tempearture after rotor
    w0 : relative velocity before rotor
    ct1 : tangential absolute velocity after rotor
    massflow : mass flow rate of air
    flowarea1 : flow area at after rotor
    spool_speed1 : linear velocity at pitchline before rotor
    spool_speed2 : linear velocity at pitchline after rotor
    cp_gas : specific heat at constant pressure
    gminusg : (gamma-1)/gamma
    R_gas : gas specific gas constant
    
    Assumptions
    -----------
    0:
    """
    cz1 = np.sqrt(velc1**2-velct1**2)
    w1 = np.sqrt(cz1**2 + (velct1-spool_speed2)**2)
    Ttr0 = Ts0  + w0**2/2/cp_gas
    Ttr1 = Ttr0 + (spool_speed2**2-spool_speed1**2)/2/cp_gas
    Ts1  = Ttr1 - w1**2/2/cp_gas
    Tt1  = Ts1  + velc1**2/2/cp_gas
    Ps1 = Pt1*(Ts1/Tt1)**(1/gminusg)
    densitys1 = Ps1/R_gas/Ts1
    cz12 = massflow/densitys1/flowarea1

    return (cz1 - cz12, Tt1, densitys1, Ts1)

def _sonic_area_ratio(M, gamma=1.4):
    '''Returns A*/A from incident Mach number'''
    top = (gamma+1)/2
    bottom = 1+(gamma-1)/2*M**2
    AstarA = M*(top/bottom)**((gamma+1)/(2*(gamma-1)))
    return AstarA

def inlet_design(Ts0, Ps0, mach0, massflow0, A0AHL, mach_throat, mach1, IPR, 
                 diffuser_angle, gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Calculates the basic parameters of the inlet based on the size of
    the incoming stream tube and certain use designtaed design
    parameters. Station numbering coincides with overall engine
    numbering
    
    Returns
    -------
    Temperature pressures and mach numbers are length 3 arrays that list
    station 1, 2 values 0: diameters of streamtube highlight, throat
    diameter,
    and fan in an array
    1: diffuser length
    2: diffuser length to throat diameter ratio
    3: static temperatures
    4: total temperatures
    5: static pressures
    6: total pressures
    7: mach numbers

    Parameters
    ----------
    Ts0 : static temperature of the freestream
    Ps0 : static pressure of the freestream
    massflow0 : mass flow rate going into the inlet
    A0AHL : ratio of the area of the freestream tube to the area of the
    inlet highlight
    mach_throat : mach number at the throat of the inlet
    mach1 : mach number at the compressor or fan face
    diffuser_angle : full angle of the diffuser that starts at the inlet
    throat and ends at the compressor or fan face

    Assumptions
    -----------
    0: No upwash induced by wings
    1: The diffuser is a frustum of a cone
    """
    
    # Freestream
    densitys0 = Ps0/R_gas/Ts0
    Tt0 = isenf.total_temperature(mach0, Ts0)
    Pt0 = isenf.total_pressure(mach0, Ps0)
    velocity0 = mach0*np.sqrt(gamma*R_gas*Ts0)
    streamtube_area = massflow0/(velocity0*densitys0)
    streamtube_diameter = np.sqrt(4*streamtube_area/np.pi)

    # Higlight
    highlight_area = streamtube_area/A0AHL
    highlight_diameter = np.sqrt(4*highlight_area/np.pi)

    # Throat
    throat_static_temp  = isenf.static_temperature(mach_throat, Tt0)
    throat_static_press = isenf.static_pressure(mach_throat, Pt0)
    speed_of_sound = np.sqrt(gamma*R_gas*throat_static_temp)
    throat_velocity = speed_of_sound*mach_throat
    throat_density = throat_static_press/(R_gas*throat_static_temp)
    throat_area = massflow0/(throat_velocity*throat_density)
    throat_diameter = np.sqrt(4*throat_area/np.pi)

    # Compressor or fan face
    fan_static_temp  = isenf.static_temperature(mach1, Tt0)
    fan_static_press = isenf.static_pressure(mach1, Pt0)
    speed_of_sound = np.sqrt(gamma*R_gas*fan_static_temp)
    fan_velocity = speed_of_sound*mach1
    fan_density = fan_static_press/(R_gas*fan_static_temp)
    fan_area = massflow0/(fan_velocity*fan_density)
    fan_diameter = np.sqrt(4*fan_area/np.pi)

    # Diffuser
    throat_radius = throat_diameter/2
    fan_radius = fan_diameter/2
    diffuser_length = (fan_radius-throat_radius)/(
        np.tan(diffuser_angle*np.pi/180))
    
    Pt1 = Pt0*IPR
    Ps1 = isenf.static_pressure(mach1, Pt1)
    Tt1 = Tt0
    Ts1 = isenf.static_temperature(mach1, Tt1)

    static_temperatures = np.array([Ts0, Ts1])
    total_temperatures = np.array([Tt0, Tt1])
    static_pressures = np.array([Ps0, Ps1])
    total_pressures = np.array([Pt0, Pt1])
    mach_numbers = np.array([mach0, mach1])
    diameters = np.array([streamtube_diameter, highlight_diameter, 
                          throat_diameter, fan_diameter])

    return (diameters, diffuser_length, diffuser_length/throat_diameter, 
            static_temperatures, total_temperatures, static_pressures, 
            total_pressures, mach_numbers)

def nacelle_design(mach0, mach1, area0area1, pressure_coeff, gamma=1.4):
    """
    Notes
    -----
    This function calculates the area ratio of external cowl area that
    produces sonic flow on the nacelle to highlight area, if the
    critical pressure coefficient is given. The flow that does not enter
    the inlet accelerates as it flows around the outside of the cowl and
    if the cowl is too large the flow will become supersonic potentially
    inducing shock and increasing drag, similar to divergence drag on an
    airfoil. Station numbering 0: freestream, 1: highlight p375
    
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

def supersonic_pressure_recovery(mach0, std='MIL'):
    """
    Notes
    -----
    Pressure recovery of supersonic inlets according to two empirical
    models 1: military standard Mil-E_5008B and 2: AIA (Aircraft
    industries Association). Both methods are for above mach 1. Mil
    standard is piecewise mach 1-5 and 5+ Both standards are considered
    conservative (worst case pressure recovery) by todays standards
    p402
    
    Returns
    -------
    ratio of total freestream pressure to fan or compressor face total
    pressure
    
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

def compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, 
                      total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, 
                      Tt31, Pt31, massflow31, mach31, 
                      gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Calculates average blade properties, i.e. gap between blades blade
    height, across entire compressor and basic compressor properties
    Only calculates aveaged properties at pitchline of blades
    Station numbering coincides with overall engine numbering
    
    Returns
    -------
    Temperature pressures and mach numbers are length 3 arrays that list
    station 2, 3, 3.1 values
    0: angular velocity of spool
    1: number of stages
    2: length of compressor
    3: array of avgerage gap, avgerage blade height, avgerage blade
    width, inlet hub diameter, outlet hub diameter
    4: static temperatures
    5: total temperatures
    6: static pressures
    7: total pressures
    8: mach numbers
    
    Parameters
    ----------
    max_tip_diam : maximum tip diameter of the blades
    max_tip_speed : maximum tip speed of rotors 
    aspect_ratio : ratio of height to width of the blades
    work_coeff : asd
    total_work : total amount of work required be performed by the
    compressor
    inlet_radius_ratio : ratio  of blade hub radius to blade tip radius
    Tt2 : total temperature at face of compressor
    Pt2 : total pressure at face of compressor
    massflow2 : air mass flow into compressor
    Tt31 : total temperature at exit of compressor
    Pt31 : total temperature at exit of compressor
    massflow31 : air mass flow at the exit of compressor after bleed air
    is removed
    mach31 : mach number at exit of compressor

    Assumptions
    -----------
    0: inner compressor diameter is constant
    1: aspect ratio is constant for all blades
    2: pitchline of blade is linear across the compressor
    """
    # total cross sectional area of the compressor
    total_area = np.pi*max_tip_diam**2/4 
    spool_speed = (max_tip_speed)/(max_tip_diam/2) # rad/s
    spool_speed_rpm = spool_speed*(1/(2*np.pi))*(60/1) # rpm
    gap_width_ratio = .25 # gap width to blade width ratio

    densityt2 = Pt2/(R_gas*Tt2)

    Ts31 = isenf.static_temperature(mach31, Tt31)
    velocity31 = mach31*np.sqrt(gamma*R_gas*Ts31)
    densityt31 = Pt31/(R_gas*Tt31)
    densitys31 = isenf.static_density(mach31, densityt31)

    inlet_blade_height = max_tip_diam/2*(1-inlet_radius_ratio)
    inlet_hub_radius = max_tip_diam/2*inlet_radius_ratio
    inlet_hub_area = np.pi*inlet_hub_radius**2
    inlet_flow_area = total_area - inlet_hub_area

    outlet_flow_area = massflow31/(densitys31*velocity31) # area31
    outlet_hub_area = total_area - outlet_flow_area
    outlet_hub_diam = np.sqrt(outlet_hub_area*4/np.pi)
    outlet_blade_height = (max_tip_diam - outlet_hub_diam)/2

    avg_blade_height = (inlet_blade_height + outlet_blade_height)/2
    avg_blade_width = avg_blade_height/aspect_ratio
    avg_gap = gap_width_ratio*avg_blade_width
    avg_pitch_diam2  = (max_tip_diam + inlet_hub_radius*2)/2
    avg_pitch_diam31 = (max_tip_diam + outlet_hub_diam)/2
    avg_pitch_diam = (avg_pitch_diam2 + avg_pitch_diam31)/2
    avg_velocity = spool_speed * avg_pitch_diam/2

    stage_work = work_coeff*avg_velocity**2/2 # work per stage
    num_stages = total_work/stage_work
    num_stages = np.ceil(num_stages)
    compressor_length = 2*num_stages*avg_blade_width + (2*num_stages-1)*avg_gap 
    
    mach2 = rootfind.newton2(_find_mach, .469, massflow=massflow2, 
                             densityt=densityt2, Tt=Tt2, area=inlet_flow_area, 
                             gamma=gamma, R_gas=R_gas)
    Ts2 = isenf.static_temperature(mach2, Tt2)
    Ps2 = isenf.static_pressure(mach2, Pt2)
    densitys2 = isenf.static_density(mach2, densityt2)
    velocity2 = mach2*np.sqrt(gamma*R_gas*Ts2)
    
    Ts3 = Ts31
    Tt3 = Tt31
    Ps31 = isenf.static_pressure(mach31, Pt31)
    Ps3 = Ps31
    Pt3 = Pt31
    mach3 = mach31
    static_temperatures = np.array([Ts2, Ts3, Ts31])
    total_temperatures = np.array([Tt2, Tt3, Tt31])
    static_pressures = np.array([Ps2, Ps3, Ps31])
    total_pressures = np.array([Pt2, Pt3, Pt31])
    mach_numbers = np.array([mach2, mach3, mach31])
    blades_values = np.array([avg_gap, avg_blade_height, avg_blade_width, 
                              inlet_hub_radius*2, outlet_hub_diam])

    return (spool_speed_rpm, num_stages, compressor_length, blades_values, 
            static_temperatures, total_temperatures, static_pressures, 
            total_pressures, mach_numbers)

def compressor_blade_design(Tt0, Pt0, mach0, massflow, flowarea1, CPR, 
                            num_stages, stage_eff, loss_coeff_rotor, 
                            loss_coeff_stator, spool_speed, pitch_diam1, 
                            pitch_diam2, gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Calculations are acrosss one stage of a compressor and at the
    pitchline
    Station numbering coincides with compressor stage numbering: 0
    before rotor, 1 between rotor and stator, 2 after stator
    Velocity representation
    c : absolute velocity
    w : relative to the rotor
    w = wz + wtj = axial flow + azimuthal flow
    c = cz + ctj = axial flow + azimuthal flow
 
    Returns
    -------
    0:
    1:
    2:
    3:
    
    Parameters
    ----------
    Tt0 : total temperature after rotor
    Pt0 : total pressure after rotor
    mach0 : mach number at comprssor face
    massflow : mass flow rate
    flowarea1 : flow area after rotor
    CPR : compressor pressure ratio
    num_stages : number of stages
    stage_eff : target efficiency of the stage
    loss_coeff_rotor : cascade loss coefficient of the rotors
    loss_coeff_stator : cascade loss coefficient of the stators
    spool_speed : angular velocity of the spool
    pitch_diam1 : pitch diameter before the rotor
    pitch_diam2 : pitch diameter after the rotor

    Assumptions
    -----------
    0: does not consider under turning
    1: angle of attack before rotor and after stator are 0
    """
    gminusg = (gamma-1)/gamma
    cp_gas = gamma*R_gas/(gamma-1)

    spool_speed1 = spool_speed*pitch_diam1/2
    spool_speed2 = spool_speed*pitch_diam2/2
    # Before rotor
    Ts0 = isenf.static_temperature(mach0, Tt0)
    Ps0 = isenf.static_pressure(mach0, Pt0)
    densitys0 = Ps0/R_gas/Ts0
    sound_speed1 = np.sqrt(gamma*R_gas*Ts0)
    velc0 = mach0*sound_speed1
    wvel0 = np.sqrt(velc0**2 + spool_speed1**2)
    # Targets for design
    stage_press_ratio = CPR**(1/num_stages)
    stage_temp_ratio = 1 + 1/stage_eff*(stage_press_ratio**gminusg-1)
    work = cp_gas*Tt0*(stage_temp_ratio-1)
    velct1 = work/spool_speed2

    Tt2 = Tt0*stage_temp_ratio
    Pt2_target = Pt0*stage_press_ratio
    Pt1 = Pt2_target
    
    # Iterating cz1 before stator and checking if total pressure after
    # stator matches 
    velc1 = 522.87*.0254*12 # initial guess
    Pt2 = 0 # initial
    while abs(Pt2_target-Pt2) > 1e-6:
        velc1 = rootfind.newton2(_find_cz, velc1, Ts0=Ts0, Pt1=Pt1, w0=wvel0, 
                                 velct1=velct1, massflow=massflow, 
                                 flowarea1=flowarea1, 
                                 spool_speed1=spool_speed1, 
                                 spool_speed2=spool_speed2, 
                                 cp_gas=cp_gas, gminusg=gminusg, R_gas=R_gas)
        
        unused, Tt1, densitys1, Ts1 = _find_cz(velc1, Ts0, Pt1, wvel0, velct1, 
                                               massflow, flowarea1, 
                                               spool_speed1, spool_speed2, 
                                               cp_gas, gminusg, R_gas)
        
        Pt1s = Pt0*(Tt1/Tt0)**(1/gminusg) # isentropic total pressure
        pressureloss_rotor = loss_coeff_rotor*(0.5*densitys0*wvel0**2)
        Pt1 = Pt1s - pressureloss_rotor
        pressureloss_stator = loss_coeff_stator*(0.5*densitys1*velc1**2)
        Pt2 = Pt1 - pressureloss_stator

        Tt2 = Tt2*(Pt2/Pt2_target)**-gminusg
        work = cp_gas*(Tt2-Tt0)
        velct1 = work/spool_speed2

    mach1 = velc1/np.sqrt(gamma*R_gas*Ts1)
    Ps1 = isenf.static_pressure(mach1, Pt1)
    densityt2 = Pt2/R_gas/Tt2
    mach2 = rootfind.newton2(_find_mach, 0.1, massflow=massflow, 
                             densityt=densityt2, Tt=Tt2, area=flowarea1, 
                             gamma=gamma, R_gas=R_gas)
    Ts2 = isenf.static_temperature(mach2, Tt2)
    Ps2 = isenf.static_pressure(mach2, Pt2)
    velc2 = mach2*np.sqrt(gamma*R_gas*Ts2)
    reaction = (Ps1 - Ps0)/(Ps2 - Ps0)

    stage_eff = ((Pt2/Pt0)**gminusg-1)/(Tt2/Tt0 - 1)

    velcz1 = np.sqrt(velc1**2-velct1**2)
    velwt1 = velct1 - spool_speed2
    
    alpha0 = 0
    alpha1 = np.arctan(velct1/velcz1)
    alpha2 = 0

    beta0 = np.arctan(0/velc0)
    beta1 = np.arctan(velwt1/velcz1)
    beta2 = alpha2

    static_temperatures = np.array([Ts0, Ts1, Ts2])
    total_temperatures  = np.array([Tt0, Tt1, Tt2])
    static_pressures    = np.array([Ps0, Ps1, Ps2])
    total_pressures     = np.array([Pt0, Pt1, Pt2])

    absolute_angles = np.array([alpha0, alpha1, alpha2])
    relative_angles = np.array([beta0, beta1, beta2])

    velw0 = velc0
    velw1 = np.sqrt(velcz1**2+velwt1**2)
    velw2 = velc2
    absolute_velocities = np.array([velc0, velc1, velc2])
    relative_velocities = np.array([velw0, velw1, velw2])
    mach_numbers        = np.array([mach0, mach1, mach2])

    return (reaction, stage_eff, static_temperatures, total_temperatures, 
            static_pressures, total_pressures, absolute_angles, 
            relative_angles, absolute_velocities, relative_velocities, 
            mach_numbers)

def airfoil_layout(beta0, beta1, alpha1, alpha2, rotor_solidity, rotor_diam, 
                   rotor_width, stator_solidity, stator_diam, stator_width, 
                   max_thick, le_thick, te_thick):
    """
    Notes
    -----
    Calculates how many airfoils are needed based on compressor design
    calculations and user defined design targets. Station numbering
    coincides with compressor stage numbering: 0 before rotor, 1 between
    rotor and stator, 2 after stator

    Returns
    -------
    0: number of airfoils
    1: rotor values - chord, spacing, camber angle, stagger angle
    2: stator values  - chord, spacing, camber angle, stagger angle
    3: rotor thicknesses - leading edge radius, trailing edge radius,
    maximum thickness
    4: stator thicknesses - leading edge radius trailing edge radius,
    maximum thickness

    Parameters
    ----------
    beta0 : angle between flow entering rotor and centerline of engine
    beta1 : angle between flow exiting  rotor and centerline of engine
    alpha1 : angle between flow entering stator and centerline of engine
    alpha2 : angle between flow entering stator and centerline of engine
    rotor_solidity : ratio of rotor chord to rotor spacing
    rotor_diam : diameter at which the rotor properties exist
    rotor_width : width of the rotor when measured parallel to the engine
    stator_solidity : ratio of stator chord to stator spacing
    stator_diam : diameter at which the stator properties exist
    stator_width : width of the stator when measured parallel to the engine
    max_thick : maximum thickness of the airfoils, given in percentage
    of chord both stator and rotor will use this value
    le_thick : thickness of the leading  edge of the airfoils given in
    percentage of chord both stator and rotor will use this value
    te_thick : thickness of the trailing edge of the airfoils given in
    percentage of chord both stator and rotor will use this value

    Assumptions
    -----------
    0: does not consider under turning leaving trailing edge
    """
    
    rotor_stagger_angle = (beta0 + beta1)/2
    rotor_chord = rotor_width/np.cos(rotor_stagger_angle)
    rotor_spacing = rotor_chord/rotor_solidity
    rotor_num_airfoils  = np.ceil(np.pi*rotor_diam/rotor_spacing)
    rotor_camber_angle = beta0 - beta1

    stator_stagger_angle = (alpha1 + alpha2)/2
    stator_chord = stator_width/np.cos(stator_stagger_angle)
    stator_spacing = stator_chord/stator_solidity
    stator_num_airfoils = np.ceil(np.pi*stator_diam/stator_spacing)
    stator_camber_angle = alpha1 - alpha2

    num_airfoils = np.array([rotor_num_airfoils, stator_num_airfoils], 
                            dtype=int)
    
    stator_le_thick = le_thick*stator_chord
    stator_te_thick = te_thick*stator_chord
    stator_thick = max_thick*stator_chord

    rotor_le_thick = le_thick*rotor_chord
    rotor_te_thick = te_thick*rotor_chord
    rotor_thick = max_thick*rotor_chord

    rotor_values = np.array([rotor_chord, rotor_spacing, rotor_camber_angle, 
                             rotor_stagger_angle])
    stator_values = np.array([stator_chord, stator_spacing, stator_camber_angle, 
                              stator_stagger_angle])
    rotor_thicks = np.array([rotor_le_thick, rotor_te_thick, rotor_thick])
    stator_thicks = np.array([stator_le_thick, stator_te_thick, stator_thick])

    return np.array([num_airfoils, rotor_values, stator_values, rotor_thicks, 
                     stator_thicks])

def combustor_design(Tt31, Pt31, airflow, ref_vel, pitch_diam, flow_split, 
                     passage_vel, min_diam_casing, max_diam_casing, 
                     max_dome_vel, comblen_domeheight, fuelflow, LHV, 
                     length_height, wall_angle, height_turbine_inlet, 
                     gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Design of an annular combustor using spacerate
    Spacerate is a concept used a GE, and is heat released per second
    per volume. If the space rate is too small there is not enough heat
    released in the combustor volume and combustion will spill over to
    the turbine.
    Iteration process
    Station numbering coincides with overall engine numbering
 
    Returns
    -------
    0: number of fuel injectors
    1: diameter of inner casing
    2: diameter of outer casing
    3: diameter of inner pass
    4: diameter of outer pass
    5: combustor length
    6: inlet prediffuser length
    7: inlet prediffuser height
    8: reference height
    9: velocity of air in dome
    10: dome height
    11: spacerate

    Parameters
    ----------
    Tt31 : total tempearture entering prediffuser
    Pt31 : total pressure entering prediffuser
    airflow : mass flow rate of the air
    ref_vel : velocity of air when fully diffused
    pitch_diam : pitch diameter of the annulus
    flow_split : fraction of the air that is to be split above, inside,
    and below the combustor
    passage_vel : desired velocity of the passages above and below the
    combustor
    min_diam_casing : innermost diameter of combustor walls from
    centerline of engine
    max_diam_casing : outermost diameter of combustor walls fome
    centerline of engine
    max_dome_vel : maximum velocity inside the dome (dome is the
    structure causing recirculation)
    comblendomeheight : ratio of combustor length to height of dome
    fuelflow : mass flow rate of the fuel
    LHV : lower heating value of the fuel
    length_height : ratio of the length to height of the pre-diffuser
    wall_angle : half angle of the pre-diffuser
    height_turbine_inlet : height of passageway from combustor to
    turbine

    Assumptions
    -----------
    0: combustion annulus is centered around the pitchline
    1: injectors are centered around the pitchline
    """
    ref_vel_max = 100
    passage_vel_max = 180
    cldh_max = 2.5
    cldh_min = 2

    # Booleans for deteriming of the values are within the limits
    casing = False
    dome = False
    spacerate = False

    while casing is False:
        diam_inner_casing, diam_outer_casing, rhot31, ref_height, ref_area = (
            _combustor_casing_calcs(Tt31, Pt31, airflow, ref_vel, pitch_diam, 
                                    R_gas))
        if diam_inner_casing < min_diam_casing or diam_outer_casing > max_diam_casing:
            ref_vel += 1
            if ref_vel > ref_vel_max:
                ref_vel = ref_vel_max
                passage_vel += 1
                if passage_vel > passage_vel_max:
                    print("""ERROR: Reference velocity and passage velocity at
                    maximum and casing geomtery exceeds limits""")
                    break
        else:
            casing = True

    while dome is False:
        (dome_height, dome_vel, diam_inner_pass, diam_outer_pass) = (
            _combustor_dome_calcs(flow_split, airflow, passage_vel, 
                                  diam_inner_casing, diam_outer_casing, rhot31,
                                  ref_height, pitch_diam))
        if dome_vel > max_dome_vel:
            ref_vel -= 1
            casing = False
        else:
            dome = True
    
    while casing is False:
        diam_inner_casing, diam_outer_casing, rhot31, ref_height, ref_area = (
            _combustor_casing_calcs(Tt31, Pt31, airflow, ref_vel, pitch_diam, 
                                    R_gas))
        if diam_inner_casing < min_diam_casing or (
            diam_outer_casing > max_diam_casing):
            passage_vel += 1
            if passage_vel > passage_vel_max:
                msg = '''ERROR: Passage velocity at maximum, reference velocity
                      at limit and casing geomtery exceeds limits.'''
                print(msg)
                break
        else:
            casing = True
    
    while spacerate is False:
        space_rate, comb_length = _space_rate_calcs(Pt31, airflow, fuelflow, 
                                                    pitch_diam, LHV, 
                                                    comblen_domeheight, 
                                                    height_turbine_inlet, 
                                                    dome_height)
        if space_rate > 10e6:
            comblen_domeheight += .02
            if comblen_domeheight > cldh_max:
                comblen_domeheight = cldh_max
                break
        elif space_rate < 8e6:
            comblen_domeheight -= .02
            if comblen_domeheight < cldh_min:
                comblen_domeheight = cldh_min
                break
        else:
            spacerate = True

    while spacerate is False:
        space_rate, comb_length = _space_rate_calcs(Pt31, airflow, fuelflow, 
                                                    pitch_diam, LHV, 
                                                    comblen_domeheight, 
                                                    height_turbine_inlet, 
                                                    dome_height)

        if space_rate > 10e6:
            # Volume is too small so the dome velocity needs to be
            # decreased which increases the dome height which increases
            # combustor volume
            dome_vel -= 1*.0254*12
            dome_area = flow_split*airflow/rhot31/dome_vel
            dome_height = dome_area/np.pi/pitch_diam
        elif space_rate < 8e6:
            # Volume is too large so the dome velocity needs to be
            # increased which changes the dome height which decreases
            # combustor volume
            dome_vel += 1
            dome_area = flow_split*airflow/rhot31/dome_vel
            dome_height = dome_area/np.pi/pitch_diam
            if dome_vel > max_dome_vel:
                print("""ERROR: Dome velocity exceeded maximum and space 
                rate is too low""")
                break
        else:
            spacerate = True

    diff_area_ratio = 1+2*length_height*np.tan(wall_angle)
    inlet_area = ref_area/diff_area_ratio
    inlet_height = inlet_area/np.pi/pitch_diam
    inlet_length = inlet_height*length_height

    circumference = np.pi*pitch_diam
    num_nozzles = np.ceil(circumference/dome_height)

    return np.array([num_nozzles, diam_inner_casing, diam_outer_casing, 
                     diam_inner_pass, diam_outer_pass, comb_length, 
                     inlet_length, inlet_height, ref_height, dome_vel, 
                     dome_height, space_rate])

def combustor_pressure_loss(Tt31:float, Pt31:float, Tt4:float, area31:float, 
                            dome_area:float, mach31:float, gamma:float, 
                            dome_vel:float=None, R_gas:float=287.05, 
                            method='book'):
    """
    Notes
    -----
    This function attempts to characterize flow properties in the
    combustor. There is an empirical relationship for total presure loss
    due to the prediffuser, heat addition in the combustor, and one for
    mach number of air to be combusted. Prediffuser pressure loss is
    from the book, heat comubstion pressure is from class, and mach
    number is from the book. The class relationship for to-be-combusted
    mach number for air is in theis function but not used due to
    uncertainty, but is used in the overall design of the combustor. The
    book air mach number is higher than the class one and produces a
    higher change in total pressure. It is recommended to use the book
    mach number.
    Station numbering "in" entrance of dome
    p517

    Returns
    -------
    0: the change in total pressure from the pre combusted air to
    combusted gas mixture (Pt_combusted - Pt_not comubsted)
    1: total pressure of air entering combustor
    2: mach number of air to be combusted

    Parameters
    ----------
    Tt31 : total temperature of air leaving compressor
    Pt31 : total pressure of air leaving compressor
    Tt4 : maximum temperature of exhuast leaving the combustor
    area31 : flow area of air leaving compressor
    dome_area : area of combustor dome
    mach31 : mach number of air leaving combustor
    gamma : ratio of specific heats of the air
    dome_vel : velocity in the dome
    R_gas : gas constant of the air
    method : book, class : choose whic formula for mach number of air in
    the combustor

    Assumptions 
    -----------
    0: prediffuser design in a dump diffuser
    1: an area ratio 1<A2/A1<5 
    """

    first = (1-area31/dome_area)
    Ptin = Pt31*np.exp(-gamma*mach31**2/2*(first**2+first**6))
    
    match method:
        case 'book':
            # derived from book p517
            machin = rootfind.newton2(prediff_mach, .5, mach1=mach31, 
                                      Pt31=Pt31, Ptin=Ptin, area31=area31, 
                                      dome_area=dome_area, gamma=gamma) 
        case 'class':
            machin = dome_vel/np.sqrt(gamma*R_gas*Tt31)
        case _:
            machin = None
            print("Specifiy 'book' or 'class'")

    deltaPt = Ptin*0.53*machin**2*(0.95+.05*(Tt4/Tt31))

    return deltaPt, Ptin, machin

def turbine_blade_design(Tt0, Pt0, Tt2, Pt2, work, alpha2, massflow, tip_diam, 
                         hub_diam, spool_speed, stage_eff, gamma_hot, 
                         gamma_cold=1.4, R_gas=287.05):
    """
    Notes
    -----
    Calculations are acrosss one stage of a trubine and at the pitchline
    Station numbering coincides with compressor stage numbering: 0
    before stator, 1 between stator and rotor, 3 after rotor
    Velocity representation
    velc : absolute velocity
    velw : relative to the rotor
    velw = velwz + velwtj = axial flow + azimuthal flow
    velc = velcz + velctj = axial flow + azimuthal flow
 
    Returns
    -------
    Eveything but reaction and flow_area is a length 3 numpy array with
    the values for the stage in order
    0: reaction
    1: flow area
    2: static temperatures 
    3: total temperatures 
    4: static pressures 
    5: total pressures 
    6: absolute angles
    7: relative angles
    8: absolute velocities
    9: relative velocities
    10: mach numbers

    Parameters
    ----------
    Tt0 : total temperature in front of the stator
    Pt0 : total pressure in front of the stator
    Tt2 : total temperature after the rotor
    Pt2 : total pressure after the rotor
    work : work required by turbine
    alpha2 : angle of attack of exhaust after the turbine
    massflow : mass flow rate
    tip_diam : turbine blade tip diameter
    hub_diam : turbine blade hub diameter
    spool_speed : angular velcoity of the spool
    stage_eff : efficiency of the stage
    gamma_hot : ratio of specific heat for combusted gases

    Assumptions
    -----------
    0: single stage
    1: constant pitchline
    2: constant flow area
    3: cooling flows not considered
    4: adiabatic
    5: no cascade losses
    """
    
    pitch_diam1 = (tip_diam + hub_diam)/2
    spool_speed1 = spool_speed*pitch_diam1/2
    # Before stator
    flow_area = np.pi*(tip_diam**2 - hub_diam**2)/4
    densityt0 = Pt0/Tt0/R_gas
    mach0 = rootfind.newton2(_find_mach, .3, massflow=massflow, 
                             densityt=densityt0, Tt=Tt0, area=flow_area, 
                             gamma=gamma_hot, R_gas=R_gas)
    Ps0 = isenf.static_pressure(mach0, Pt0, gamma=gamma_hot)
    Ts0 = isenf.static_temperature(mach0, Tt0, gamma=gamma_hot)
    velc0 = mach0*np.sqrt(gamma_hot*R_gas*Ts0)

    # After stator
    velct1 = work/spool_speed1
    velwt1 = velct1 - spool_speed1

    mach1 = rootfind.newton2(_find_mach2, .9, massflow=massflow, 
                             densityt=densityt0, Tt=Tt0, flow_area=flow_area, 
                             velocity_comp=velct1, gamma=gamma_hot, 
                             R_gas=R_gas)
    unused, velcz1 = _find_mach2(mach1, massflow, densityt0, Tt0, flow_area, 
                                 velct1, gamma=gamma_hot, R_gas=R_gas)
    alpha1 = np.arctan(velct1/velcz1)
    Ps1 = isenf.static_pressure(mach1, Pt0, gamma_hot)
    Ts1 = isenf.static_temperature(mach1, Tt0, gamma_hot)
    velc1 = np.sqrt(gamma_hot*R_gas*Ts1)*mach1
    velcz1 = np.sqrt(velc1**2-velct1**2)
    beta1 = np.arctan(velwt1/velcz1)

    # After rotor
    densityt2 = Pt2/Tt2/R_gas
    mach2 = rootfind.newton2(_find_mach, .3, massflow=massflow, 
                             densityt=densityt2, Tt=Tt2, area=flow_area, 
                             gamma=gamma_hot, R_gas=R_gas)
    Ps2 = isenf.static_pressure(mach2, Pt2, gamma_hot)
    Ts2 = isenf.static_temperature(mach2, Tt2, gamma_hot)
    reaction = (Ts1-Ts2)/(Tt0-Ts2)
    velcz2 = np.sqrt(gamma_hot*R_gas*Ts2)*mach2
    velwt2 = spool_speed1 # because flow is purely axial exiting turbine
    beta2 = np.arctan(velwt2/velcz2)
    
    Tt1 = Tt0
    static_temperatures = np.array([Ts0, Ts1, Ts2])
    total_temperatures = np.array([Tt0, Tt1, Tt2])
    static_pressures = np.array([Ps0, Ps1, Ps2])
    total_pressures = np.array([Pt0, Pt2, Pt2])

    alpha0 = 0
    alpha2 = 0
    beta0 = 0
    absolute_angles = np.array([alpha0, alpha1, alpha2])
    relative_angles = np.array([beta0, beta1, beta2])

    velc2 = velcz2
    velw0 = velc2
    velw1 = np.sqrt(velcz1**2+velwt1**2)
    velw2 = np.sqrt(velcz2**2+velwt2**2)
    absolute_velocities = np.array([velc0, velc1, velc2])
    relative_velocities = np.array([velw0, velw1, velw2])
    mach_numbers = np.array([mach0, mach1, mach2])

    return (reaction, flow_area, static_temperatures, total_temperatures, 
            static_pressures, total_pressures, absolute_angles, 
            relative_angles, absolute_velocities, relative_velocities, 
            mach_numbers)

def afterburner_pressure_loss(machi, drag_coeff, gammai, gammae, q=None, 
                              cpi=None, Ti=None, mode='dry'):
    """
    Notes
    -----
    Total pressure loss in a constant area afterburner in either dry or
    wet afterburner drag coefficient notes. Subscripts i and e denote
    properties at entry and exit of afterburner. Derivation of equation
    7.97 for mache^2 is wrong and should be mache^2 =
    (A^2ge-2(ge-1)-A\sqrt((A^2+2)ge^2+2) )/(2-A^2)/(ge-1)/ge where g is
    gamma, all other equations are correct for the derivation of
    Pte/Pti. Figure 7.32 of mache vs machi reflects the incorrect book
    equation, figure 7.31 cannot be replicated, figure 7.34 cannot be
    replicated and contradicts figure 7.31. The code reflects the
    correct equation for mach2 presented here
    p506

    Returns
    -------
    ratio of total pressure at the beginning of the afterburner to the
    total pressure at the exit of the afterburner

    Parameters
    ----------
    machi : mach number at beginning of afterburner
    drag_coeff : drag coefficient of the flame holders
    gammai : ratio of specific heats at beginning of afterburner
    gammae : ratio of specific heats at exit of afterburner
    dry : no fuel combusted in afterburner
    wet : fuel combusted in afterburner
    q : amount of heat released when fuel is combusted
    cpi : specific heat at constant pressure of exhuast entering
    afterburner
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
    mache_sq = ( third**2*gammae-2*geminus-third*
                np.sqrt((third**2-2)*gammae**2+2) )/second

    sixth = 1+gammai*machi**2*(1-drag_coeff/2)
    seventh = 1+gammae*mache_sq
    pressure_ratio = sixth/seventh*(1+geminus/2*mache_sq)**(1/geminusg)/(
        1+giminus/2*machi**2)**(1/giminusg)

    return pressure_ratio

def nozzle_design(Tt5:float, Pt5:float, flowarea5:float, Ps9:float, 
                  mach5:float, Cv:float, gamma_hot:float, 
                  gamma_cold:float=1.4, R_gas:float=287.05):
    """
    Notes
    -----
    Calculates nozzle throat diameter, exit diameter, exit mach number
 
    Returns
    -------
    Temperatures and pressures are length 2 arrays for the throat and
    exit properties
    0: throat diameter
    1: exit diameter
    2: exit mach number
    3: static temperatures
    3: total temperatures
    3: static pressures
    3: total pressures
    
    Parameters
    ----------

    Assumptions
    -----------
    0: perfectly matched pressure conditions
    1: does not correct geometrical throat area to effective throat area
    """ 
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)
    gminus = gamma_hot - 1
    gminusg = gminus/gamma_hot

    AstarA = _sonic_area_ratio(mach5, gamma_hot)
    Astar = flowarea5*AstarA

    Ts9 = Tt5*(1-Cv**2*(1-(Ps9/Pt5)**gminusg))
    velocity9 = np.sqrt(2*cp_hot*(Tt5-Ts9))
    mach9 = velocity9/np.sqrt(gamma_hot*R_gas*Ts9)
    AstarAexit = _sonic_area_ratio(mach9, gamma_hot)
    Aexit = AstarAexit**-1*Astar

    throat_diameter = np.sqrt(4/np.pi*Astar)
    exit_diameter = np.sqrt(4/np.pi*Aexit)

    Ps8 = isenf.static_pressure(1, Pt5)
    Ts8 = isenf.static_temperature(1, Tt5)
    Tt8 = Tt5
    Pt8 = Pt5
    Pt9 = isenf.total_pressure(mach9, Pt5)
    Tt9 = isenf.total_temperature(mach9, Tt5)

    static_temperatures = np.array([Ts8, Ts9])
    total_temperatures = np.array([Tt8, Tt9])
    static_pressures = np.array([Ps8, Ps9])
    total_pressures = np.array([Pt8, Pt9])

    return (throat_diameter, exit_diameter, mach9, static_temperatures, 
            total_temperatures, static_pressures, total_pressures)

def inlet_design_test_values():
    Ps0 = 4.36 # psia
    Ts0 = 411.5 # R
    mach0 = .8
    massflow2 = 70.1 #lbm/s
    A0Ahl = .7
    throat_mach = .7 # max is 0.75
    fan_mach = .4
    diffuser_angle = 5
    IPR = .99
    massflow2 = units.convert_mass(massflow2)
    Ps0 = Ps0*6895
    Ts0 = Ts0*5/9

    inlet_values = inlet_design(Ts0, Ps0, mach0, massflow2, A0Ahl, throat_mach,
                                 fan_mach, IPR, diffuser_angle)

def compressor_design_test_values():
    """Compressor Design"""
    Tt2 = 464.5 # R
    Pt2 = 6.58 # psia
    massflow2 = 70.17 # lbm/s
    Tt31 = 897.7 # R
    Pt31 = 52.67 # psia
    mach31 = .3
    massflow31 = 63.15 # lbm/s
    work_coeff = .6
    total_work = 103.97 # Btu/(lbm/sec)
    comp_press_ratio = 8
    comp_eff = .87
    aspect_ratio = 2.5 # height to width ratio
    inlet_radius_ratio = .4 # hub radius to tip radius ratio
    max_tip_diam = 29.58 # in
    max_tip_speed = 1500 # ft/s

    Tt2, Tt31 = units.convert_temperature([Tt2, Tt31], 'K')
    Pt2, Pt31 = units.convert_pressure([Pt2, Pt31], 'Pa')
    massflow2 = massflow2/2.205
    massflow31 = massflow31/2.205
    max_tip_diam *= .0254 # m
    max_tip_speed *= 12*.0254 # m/s
    total_work = units.convert_energy(total_work, 'J') # J/(kg/s)

    compressor_values = compressor_design(max_tip_diam, max_tip_speed, 
                                          aspect_ratio, work_coeff, total_work,
                                            inlet_radius_ratio, Tt2, Pt2, 
                                            massflow2, Tt31, Pt31, massflow31, 
                                            mach31)

def compressor_blade_test_values():
    Tt1 = 464.5 * 5/9
    Pt1 = 6.58 *6895
    massflow = 70.17 / 2.20462
    flowarea2 = 551.19 *.0254**2 # in^2
    loss_coeff_rotor = .1
    loss_coeff_stator = .07
    CPR = 8
    num_stages = 6
    stage_eff = 0.87
    mach0 = .50245
    spool_speed = 11619*2*np.pi/60
    pitch_diam1 = 20.71*.0254
    pitch_diam2 = 21.37*.0254

    cmp_bld_vals = compressor_blade_design(Tt1, Pt1, mach0, massflow, 
                                           flowarea2, CPR, num_stages, 
                                           stage_eff, loss_coeff_rotor, 
                                           loss_coeff_stator, spool_speed, 
                                           pitch_diam1, pitch_diam2)

    rotor_solidity = 1.3 # chord/spacing
    rotor_diam = 20.75 * .0254
    rotor_width = 3.56 *.0254

    stator_solidity = 1.2 # chord/spacing
    stator_diam = 22*.0254  
    stator_width = 3.3 *.0254

    meanline_slope = 5.9 # deg
    max_thick = .05
    le_thick = .009
    te_thick = .009

    beta0 = -59.4*np.pi/180
    beta1 = 43.3*np.pi/180
    alpha1 = 25.2*np.pi/180
    alpha2 = 0*np.pi/180

    airfoil_layout(beta0, beta1, alpha1, alpha2, rotor_solidity, rotor_diam, 
                   rotor_width, stator_solidity, stator_diam, stator_width,
                   max_thick, le_thick, te_thick)

def combustor_design_test_values():
    # Compressor parameters
    Tt31 = 897.1 # R
    Pt31 = 52.64 # psia
    airflow = 62.86 # lbm/s
    pitch_diam = 28.02 # in
    tip_diam = 29.51 # in
    # Combustor parameters
    Tt4 = 2560 # R
    Pt4 = 50.01 # psia
    fuelflow = 1.55 # lbm/s
    LHV = 18550 # Btu/lbm
    max_diam_casing = 42 # in
    min_diam_casing = 13 # in
    lengthheight = 2
    wall_angle = 5 # deg
    flow_split = 1/3
    # Turbine Parameters
    height_turbine_inlet = 4
    # Initial parameters
    ref_vel = 90 # ft/sec max = 100
    passage_vel = 150 # ft/sec max=180
    dome_vel_max = 80 # ft/sec
    comblendomeheight = 2.25 # ratio max = 2.5 min=2
    # Conversions
    Tt31, Tt4 = units.convert_temperature([Tt31, Tt4])
    Pt31, Pt4 = units.convert_pressure([Pt31, Pt4])
    airflow, fuelflow = units.convert_mass([airflow, fuelflow])
    pitch_diam *= .0254
    tip_diam *=.0254
    max_diam_casing *= .0254
    min_diam_casing *= .0254
    wall_angle *= np.pi/180
    ref_vel *= .0254*12
    dome_vel_max *= .0254*12
    passage_vel *= .0254*12
    height_turbine_inlet *= .0254

    combustor_values = combustor_design(Tt31, Pt31, airflow, ref_vel, 
                                        pitch_diam, flow_split, passage_vel, 
                                        min_diam_casing, max_diam_casing, 
                                        dome_vel_max, comblendomeheight, 
                                        fuelflow, LHV, lengthheight, 
                                        wall_angle, height_turbine_inlet)
    
    return combustor_values
   
def turbine_blade_test_values():
    Ps0 = 4.364
    Tt4 = 2560 # R
    Pt4 = 50.04 # psia
    massflow4 = 64.71 # lbm/s
    N = 11619 # rpm
    eff_turb = .92 # must stay above .92
    Tt49 = 2149 # R
    Pt49 = 23.23 # psia
    gamma_hot = 4/3
    work = 112.74
    mach49 = .55
    alpha0 = 0
    alpha2 = 0 # no exit swirl

    tip_diam = 32.55 # in
    hub_diam = 25.76 # in
    nozzle_coeff = .983
    Tt4  = units.convert_temperature(Tt4)
    Pt4  = units.convert_pressure(Pt4)
    Ps0  = units.convert_pressure(Ps0)
    Pt49 = units.convert_pressure(Pt49)
    Tt49 = units.convert_temperature(Tt49)
    massflow4 = units.convert_mass(massflow4)
    work = units.convert_energy(work)
    tip_diam = tip_diam *.0254
    hub_diam = hub_diam *.0254
    N *= 2*np.pi/60

    turbine_values = turbine_blade_design(Tt4, Pt4, Tt49, Pt49, work, alpha2, 
                                          massflow4, tip_diam, hub_diam, N, 
                                          eff_turb, gamma_hot)
    
    nozzle_values = nozzle_design(turbine_values[3][2], turbine_values[5][2], 
                                  turbine_values[1], Ps0, 
                                  turbine_values[10][2], nozzle_coeff, 
                                  gamma_hot)

if __name__ == '__main__':
    # inlet_design_test_values()
    # compressor_design_test_values()
    # compressor_blade_test_values()
    combustor_design_test_values()
    # turbine_blade_test_values()
    

