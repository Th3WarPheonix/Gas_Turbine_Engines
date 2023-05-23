
import numpy as np
import IsentropicFlow as isenf
import RootFinding as rootfind
import matplotlib.pyplot as plt

def convert_temps(temps, to:str='K'):
    """Convert temperature between Kelvin and Rankine
    to = 'K' for converting to Kelvin
    to = 'R' for converting to Rankine"""
    con_factor = 1.8
    if to == 'R':
        try:
            temps *= con_factor
            return temps
        except:
            return np.array(temps)*con_factor
    elif to == 'K':
        try:
            temps = temps/con_factor
            return temps
        except:
            return np.array(temps)/con_factor
    else:
        print('Did not convert temperatures')

def convert_mass(mass, to:str='kg'):
    """Convert mass between lbm and kg
    to = 'lbm' for converting to lbm
    to = 'kg' for converting to kg"""
    con_factor = 2.20462
    if to == 'lbm':
        try:
            mass *= con_factor
            return mass
        except:
            return np.array(mass)*con_factor
    elif to == 'kg':
        try:
            mass = mass/con_factor
            return mass
        except:
            return np.array(mass)/con_factor
    else:
        print('Did not convert masses')

def convert_pressures(pressures, to:str='Pa'):
    """Convert pressure between Pascals and PSI
    to = 'Pa' for converting to Pascals
    to = 'psi' for converting to PSI"""
    con_factor = 6895.0
    if to == 'psi':
        try:
            pressures = pressures/con_factor
            return pressures
        except:
            return np.array(pressures)/con_factor
    elif to == 'Pa':
        try:
            pressures = pressures*con_factor
            return pressures
        except:
            return np.array(pressures)*con_factor
    else:
        print('Did not convert pressures')

def convert_energy(works, to:str='J'):
    """Convert mass specific energy/work between Btu/lbm and J/kg
    to = 'J' for converting to J/kg
    to = 'BTU' for converting to BTU/lbm"""
    if to == 'BTU':
        try:
            works = works / 1055 / 2.205
            return works
        except:
            return np.array(works) / 1055 / 2.205
    elif to == 'J':
        try:
            works = works * 1055 * 2.205
            return works
        except:
            return np.array(works) * 1055 * 2.205
    else:
        print('Did not convert energy')
        
def inlet_design(stream_density:float, stream_velocity:float, massflow:float, A0AHL:float, mach_throat:float, Tt0:float, Pt0:float, mach2:float, diffuser_angle:float, gamma:float=1.4, R_gas:float=287.05):
    """
    Notes
    -----
    Calculates the basic parameters of the inlet based on the size of the incoming stream tube and certain use designtaed design parameters
    Station numbering coincides with overall engine numbering
    
    Returns
    -------
    0: stream tube diameter, highlight diameter, throat diameter, fan diameter
    1: diffuser length
    2: diffuser length to throat diameter ratio
    
    Parameters
    ----------
    stream_density  : density of the freestream
    stream_velocity : velocity of the freestream relative to the inlet
    massflow : asd
    A0AHL : ratio of the area of the freestream tube to the area of the inlet highlight
    mach_throat : mach number at the throat of the inlet
    Tt0 : total temperature of the freestream
    Pt0 : total pressure of the freestream
    Ts1 : asd
    Ps1 : asd
    mach2 : mach number at the compressor or fan face
    diffuser_angle : full angle of the diffuser that starts at the inlet throat and ends at the compressor or fan face

    Assumptions
    -----------
    0: No upwash induced by wings
    1: The diffuser is a frustum of a cone
    """
    
    '''Entrance calcs'''
    streamtube_area = massflow/(stream_velocity*stream_density)
    streamtube_diameter = np.sqrt(4*streamtube_area/np.pi)
    highlight_area = streamtube_area/A0AHL
    highlight_diameter = np.sqrt(4*highlight_area/np.pi)

    '''Throat calcs'''
    throat_static_temp  = isenf.T(mach_throat, Tt0)
    throat_static_press = isenf.p(mach_throat, Pt0)
    speed_of_sound = np.sqrt(gamma*R_gas*throat_static_temp)
    throat_velocity = speed_of_sound*mach_throat
    throat_density = throat_static_press/(R_gas*throat_static_temp)
    throat_area = massflow/(throat_velocity*throat_density)
    throat_diameter = np.sqrt(4*throat_area/np.pi)

    '''Compressor face or fan face calcs'''
    fan_static_temp  = isenf.T(mach2, Tt0)
    fan_static_press = isenf.p(mach2, Pt0)
    speed_of_sound = np.sqrt(gamma*R_gas*fan_static_temp)
    fan_velocity = speed_of_sound*mach2
    fan_density = fan_static_press/(R_gas*fan_static_temp)
    fan_area = massflow/(fan_velocity*fan_density)
    fan_diameter = np.sqrt(4*fan_area/np.pi)

    '''Diffuser length calcs'''
    throat_radius = throat_diameter/2
    fan_radius = fan_diameter/2
    diffuser_length = (fan_radius-throat_radius)/np.tan(diffuser_angle*np.pi/180)

    return np.array([streamtube_diameter, highlight_diameter, throat_diameter, fan_diameter]), diffuser_length, diffuser_length/throat_diameter

def find_mach(mach, massflow, densityt, Tt, area, gamma=1.4, R=287.05):
    """Find Mach number from given massflow, total temperature, density, and flow area"""
    densitys = isenf.static_density(mach, densityt, gamma=gamma)
    Ts = isenf.static_temperature(mach, Tt, gamma=gamma)
    sound_speed = np.sqrt(gamma*R*Ts)
    return (massflow/densitys/mach/sound_speed - area,)

def find_mach2(mach, massflow, densityt, Tt, flow_area, velocity_comp, gamma=1.4, R_gas=287.05):
    """Find Mach number from given massflow, total temperature, density, and flow area and a given velocity component"""
    Ts1 = isenf.static_temperature(mach, Tt, gamma=gamma)
    densitys1 = isenf.static_density(mach, densityt, gamma=gamma)
    velocity1zsq = (mach*np.sqrt(gamma*R_gas*Ts1))**2 - velocity_comp**2
    veclocity1zsq2 = (massflow/densitys1/flow_area)**2
    return (veclocity1zsq2 - velocity1zsq, np.sqrt(velocity1zsq))

def find_mach3(mach, Tt, velocity, gamma, R_gas):
    """Find mach number from given total temperature and velocity"""
    mach1 = velocity/np.sqrt(gamma*R_gas*isenf.static_temperature(mach, Tt))
    return (mach1 - mach,)

def compressor_design(max_tip_diam:float, max_tip_speed:float, aspect_ratio:float, work_coeff:float, total_work:float, inlet_radius_ratio:float, Tt2:float, Pt2:float, massflow2:float, Tt31:float, Pt31:float, massflow31:float, mach31:float, gamma:float=1.4, R_gas:float=287.05):
    """
    Notes
    -----
    Calculates average blade properties, i.e. gap between blades blade height, across entire compressor and basic compressor properties
    Only calculates aveaged properties at pitchline of blades
    Station numbering coincides with overall engine numbering
    
    Returns
    -------
    0: inlet hub diameter, outlet hub diameter, average gap between blades, average blade height
    1: spool_speed_rpm
    2: num_stages
    3: compressor_length
    
    Parameters
    ----------
    max_tip_diam : maximum tip diameter of the blades
    max_tip_speed : maximum tip speed of rotors 
    aspect_ratio : ratio of height to width of the blades
    work_coeff : asd
    total_work : total amount of work required be performed by the compressor
    inlet_radius_ratio : ratio  of blade hub radius to blade tip radius
    Tt2 : total temperature at face of compressor
    Pt2 : total pressure at face of compressor
    massflow2 : air mass flow into compressor
    Tt31 : total temperature at exit of compressor
    Pt31 : total temperature at exit of compressor
    massflow31 : air mass flow at the exit of compressor after bleed air is removed
    mach31 : mach number at exit of compressor

    Assumptions
    -----------
    0: inner compressor diameter is constant
    1: aspect ratio is constant for all blades
    2: pitchline of blade is linear across the compressor
    """
    total_area = np.pi*max_tip_diam**2/4 # total cross sectional area of the compressor
    spool_speed = (max_tip_speed)/(max_tip_diam/2) # rad/s
    spool_speed_rpm = spool_speed*(1/(2*np.pi))*(60/1) # rpm
    gap_width_ratio = .25 # gap width to blade width ratio

    densityt2 = Pt2/(R_gas*Tt2)

    Ts31 = isenf.T(mach31, Tt31)

    velocity31 = mach31*np.sqrt(gamma*R_gas*Ts31)
    densityt31 = Pt31/(R_gas*Tt31)
    densitys31 = isenf.density(mach31, densityt31)

    inlet_blade_height = max_tip_diam/2*(1-inlet_radius_ratio)
    inlet_hub_radius = max_tip_diam/2*inlet_radius_ratio 
    inlet_hub_area = np.pi*inlet_hub_radius**2
    inlet_flow_area = total_area - inlet_hub_area

    outlet_flow_area = massflow31/(densitys31*velocity31)

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
    compressor_length = 2*num_stages*avg_blade_width + (2*num_stages-1)*avg_gap # 6 rotors, 6 stators, 11 gaps in between all of the rotors and stators

    mach2 = rootfind.newton2(find_mach, .469, massflow=massflow2, densityt=densityt2, Tt=Tt2, area=inlet_flow_area)
    Ts2 = isenf.T(mach2, Tt2)
    Ps2 = isenf.p(mach2, Pt2)
    densitys2 = isenf.density(mach2, densityt2)
    velocity2 = mach2*np.sqrt(gamma*R_gas*Ts2)
    
    return ((inlet_hub_radius*2, outlet_hub_diam, avg_gap, avg_blade_height), spool_speed_rpm, num_stages, compressor_length)

def find_cz(c2, Ts1, Pt2, w1, ct2, massflow, flowarea2, spool_speed1, spool_speed2, cp_gas, gminusg, R_gas):
    cz2 = np.sqrt(c2**2-ct2**2)
    w2 = np.sqrt(cz2**2 + (ct2-spool_speed2)**2)
    Ttr1 = Ts1  + w1**2/2/cp_gas
    Ttr2 = Ttr1 + (spool_speed2**2 -spool_speed1**2)/2/cp_gas
    Ts2  = Ttr2 - w2**2/2/cp_gas
    Tt2  = Ts2  + c2**2/2/cp_gas
    Ps2 = Pt2*(Ts2/Tt2)**(1/gminusg)
    densitys2 = Ps2/R_gas/Ts2
    cz22 = massflow/densitys2/flowarea2

    return (cz2 - cz22, Tt2, densitys2)

def compressor_blade_design(Tt1, Pt1, mach1, massflow, flowarea2, CPR, num_stages, stage_eff, loss_coeff_rotor, loss_coeff_stator, spool_speed, pitch_diam1, pitch_diam2, gamma=1.4, R_gas=287.05):
    """
    Notes
    -----
    Calculations are acrosss one stage of a compressor
    Station numbering coincides with compressor stage numbering: 1 before rotor, 2 between rotor and stator, 3 after stator
    Velocity representation
    c : absolute velocity
    w : realtive to the rotor
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
    Tt1 : total temperature in front of the rotor
    Pt1 : total pressure in front of the rotor
    massflow : mass flow into the stage
    alpha1 : angle of attack of air with respect to the rotor
    press_ratio : pressure ratio of the stage
    num_stages : number of stages
    Dt1 :
    Dp1 :
    Dp2 :
    flow_area3 :
    spool_speed :
    stage_eff :
    loss_coeffr :
    loss_coeffs :
    reaction :
    alpha3 :

    Assumptions
    -----------
    0:
    1:
    2:
    """
    gminusg = (gamma-1)/gamma
    cp_gas = gamma*R_gas/(gamma-1)

    spool_speed1 = spool_speed*pitch_diam1/2
    spool_speed2 = spool_speed*pitch_diam2/2

    Ts1 = isenf.static_temperature(mach1, Tt1)
    Ps1 = isenf.static_pressure(mach1, Pt1)
    densitys1 = Ps1/R_gas/Ts1
    sound_speed1 = np.sqrt(gamma*R_gas*Ts1)
    c1 = mach1*sound_speed1
    w1 = np.sqrt(c1**2 + spool_speed1**2)

    stage_press_ratio = CPR**(1/num_stages)
    stage_temp_ratio = 1 + 1/stage_eff*(stage_press_ratio**gminusg-1)
    Tt3 = Tt1*stage_temp_ratio
    Pt3_target = Pt1*stage_press_ratio
    Pt2 = Pt3_target
    work = cp_gas*Tt1*(stage_temp_ratio-1)
    ct2 = work/spool_speed2
    
    c2 = 524.87*.0254*12
    Pt3 = 0
    while abs(Pt3_target-Pt3) > 1e-6:
        find_cz(c2, Ts1, Pt2, w1, ct2, massflow, flowarea2, spool_speed1, spool_speed2, cp_gas, gminusg, R_gas)
        c2 = rootfind.newton2(find_cz, c2, Ts1=Ts1, Pt2=Pt2, w1=w1, ct2=ct2, massflow=massflow,flowarea2=flowarea2, spool_speed1=spool_speed1, spool_speed2=spool_speed2, cp_gas=cp_gas, gminusg=gminusg, R_gas=R_gas)
        unused, Tt2, densitys2 = find_cz(c2, Ts1, Pt2, w1, ct2, massflow, flowarea2, spool_speed1, spool_speed2, cp_gas, gminusg, R_gas)
        Pt2s = Pt1*(Tt2/Tt1)**(1/gminusg)
        pressureloss_rotor = loss_coeff_rotor*(0.5*densitys1*w1**2)
        Pt2 = Pt2s - pressureloss_rotor

        pressureloss_stator = loss_coeff_stator*(0.5*densitys2*c2**2)
        Pt3 = Pt2 - pressureloss_stator
        Tt3 = Tt3*(Pt3/Pt3_target)**-gminusg
        work = cp_gas*(Tt3-Tt1)
        ct2 = work/spool_speed2

def airfoil_count():
    rotor_solidity = 1.3 # chord/spacing
    rotor_pitch_diam = 20.75 * .0254
    rotoR_gasfoil_width = 3.56 *.0254

    stator_solidity = 1.2 # chord/spacing
    stator_pitch_diam = 22*.0254  
    statoR_gasfoil_width = 3.3 *.0254

    meanline_slope = 5.9 # deg

    beta1 = -59.4*np.pi/180
    beta2 = 43.3*np.pi/180
    alpha2 = 25.2*np.pi/180
    alpha3 = 0*np.pi/180

    rotor_stagger_angle = (beta1 + beta2)/2
    rotor_chord = rotoR_gasfoil_width/np.cos(rotor_stagger_angle)
    rotor_spacing = rotor_chord/rotor_solidity
    rotor_num_airfoils  = np.ceil(np.pi*rotor_pitch_diam/rotor_spacing)
    rotor_camber_angle = beta1 - beta2

    stator_stagger_angle = (alpha2 + alpha3)/2
    stator_chord = statoR_gasfoil_width/np.cos(stator_stagger_angle)
    stator_spacing = stator_chord/stator_solidity
    stator_num_airfoils = np.ceil(np.pi*stator_pitch_diam/stator_spacing)
    stator_camber_angle = alpha2 - alpha3

    print('alpha2, alpha3 deg \t', alpha2*180/np.pi, alpha3*180/np.pi)
    print('stator chord in \t{}'.format(stator_chord/.0254))
    print('leading edge thickness in \t{}'.format(.009*stator_chord/.0254))
    print('trailing edge thickness in \t{}'.format(.009*stator_chord/.0254))
    print('camber angle \t{}'.format(stator_camber_angle*180/np.pi))
    print('thickness in \t{}'.format(.05*stator_chord/.0254))
    num_airfoils = np.array([rotor_num_airfoils, stator_num_airfoils], dtype=int)
    print(num_airfoils)
    
    return num_airfoils


def combustor(Tt31, Pt31, airflow, ref_vel, pitch_diam, flow_split, passage_vel, min_diam_casing, max_diam_casing, max_dome_vel, comblendomeheight, fuelflow, LHV, length_height, wall_angle, height_turbine_inlet, gamma=1.4, R_gas=287.05):
    rhot31 = Pt31/R_gas/Tt31
    ref_area = airflow/rhot31/ref_vel
    ref_height = ref_area/np.pi/pitch_diam
    diam_inner_casing = (pitch_diam/2 - ref_height/2)*2
    diam_outer_casing = (pitch_diam/2 + ref_height/2)*2
    
    if diam_inner_casing < min_diam_casing or diam_outer_casing > max_diam_casing:
        print('---ERROR: CASING DIAMETER EXCEEDS LIMITS---')
        print('Recommended action: change reference velocity')

    area_passage = flow_split*airflow/rhot31/passage_vel  
    diam_inner_pass = 2*np.sqrt(diam_inner_casing**2/4 + area_passage/np.pi)
    diam_outer_pass = 2*np.sqrt(diam_outer_casing**2/4 - area_passage/np.pi)
    height_inner_pass = (diam_inner_pass - diam_inner_casing)/2
    height_outer_pass = (diam_outer_casing - diam_outer_pass)/2

    dome_height = ref_height - height_inner_pass - height_outer_pass
    dome_area = np.pi*dome_height*pitch_diam
    dome_vel = flow_split*airflow/rhot31/dome_area
    
    if dome_vel > max_dome_vel:
        print('---ERROR: DOME VELOCITY EXCEEDS MAXIMUM---')
        print('Recommended action: reduce reference velocity')
        print('Dome Vel: {} \t Max: {}'.format(round(dome_vel, 2), round(max_dome_vel, 2)))

    comb_length = comblendomeheight * dome_height
    area_entr = ref_height * comb_length/2
    area_exit = 1/2*(ref_height+height_turbine_inlet)
    comb_vol = (area_entr+area_exit)*np.pi*pitch_diam*comb_length

    fuel_air = fuelflow/airflow
    airflow_lb = convert_mass(airflow, 'lbm')
    Ps31_atm = Pt31/101300
    comb_vol_ft3 = comb_vol/(.0254*12)**3
    space_rate = 3600*fuel_air*airflow_lb*LHV/Ps31_atm/comb_vol_ft3

    if not 8e6 < space_rate < 10e6:
        print('---ERROR: SPACE RATE EXCCEEDS LIMITS---')
        print('Recommended action: change length to height ratio, then dome velocity')
        print('Dome Vel: {:.3e} \t Bounds: {:.3e},{:.3e}'.format(space_rate, 8e6, 10e6))

    diff_area_ratio = 1+2*length_height*np.tan(wall_angle)
    inlet_area = ref_area/diff_area_ratio
    inlet_height = inlet_area/np.pi/pitch_diam
    inlet_length = inlet_height*length_height

    circumference = np.pi*pitch_diam
    num_nozzles = np.ceil(circumference/dome_height)

    return np.array((num_nozzles, diam_inner_casing, diam_outer_casing, diam_inner_pass, diam_outer_pass, comb_length, inlet_length, inlet_height, ref_height, dome_height, height_inner_pass, height_outer_pass))


def assignment5(dfConfigs, R_gas=287.05):
    Tt0 = convert_temps(dfConfigs['Config 1'].loc['0 Freestream', 'Total Temperature (R)'], 'SI')
    Pt0 = convert_pressures(dfConfigs['Config 1'].loc['0 Freestream', 'Total Pressure (psia)'], 'SI')
    Ts0 = convert_temps(dfConfigs['Config 1'].loc['0 Freestream', 'Static Temperature (R)'], 'SI')
    Ps0 = convert_pressures(dfConfigs['Config 1'].loc['0 Freestream', 'Static Pressure (psia)'], 'SI')
    Ts1 = convert_temps(dfConfigs['Config 1'].loc['1 Fan Inlet', 'Static Temperature (R)'], 'SI')
    Ps1 = convert_pressures(dfConfigs['Config 1'].loc['1 Fan Inlet', 'Static Pressure (psia)'], 'SI')

    '''Inlet Design'''
    density0 = Ps0/(R_gas*Ts0)
    velocity0 = 795.5 # ft/s
    velocity0 = velocity0 * 12*.0254
    massflow0 = 70.17 # lbm/s
    massflow0 = massflow0/2.205
    A0Ahl = .7
    throat_mach = .7
    fan_mach = .4
    diffuser_angle = 5
    diams, diff_length, diff_LH = inlet_design(density0, velocity0, massflow0, A0Ahl, throat_mach, Tt0, Pt0, fan_mach, diffuser_angle)
    return (diams/.0254, diff_length/.0254, diff_LH)

def assignment6():
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

    Tt2, Tt31 = convert_temps([Tt2, Tt31], 'SI')
    Pt2, Pt31 = convert_pressures([Pt2, Pt31], 'SI')
    massflow2 = massflow2/2.205
    massflow31 = massflow31/2.205
    max_tip_diam *= .0254 # m
    max_tip_speed *= 12*.0254 # m/s
    total_work = convert_energy(total_work, 'SI') # J/(kg/s)

    return compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, Tt31, Pt31, massflow31, mach31)

def assignment7():
    Tt1 = 464.5 * 5/9
    Pt1 = 6.58 *6895
    massflow = 70.17 / 2.20462
    flowarea2 = 551.19 *.0254**2 # in^2
    loss_coeff_rotor = .1
    loss_coeff_stator = .07
    CPR = 8
    num_stages = 6
    stage_eff = 0.87
    mach1 = .50245
    spool_speed = 11619*2*np.pi/60
    pitch_diam1 = 20.71*.0254
    pitch_diam2 = 21.37*.0254
    compressor_blade_design(Tt1, Pt1, mach1, massflow, flowarea2, CPR, num_stages, stage_eff, loss_coeff_rotor, loss_coeff_stator, spool_speed, pitch_diam1, pitch_diam2)


def assignment8():
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
    comblendomeheight = 2.25 # ratio max = 2.5
    # Conversions
    Tt31, Tt4 = convert_temps([Tt31, Tt4])
    Pt31, Pt4 = convert_pressures([Pt31, Pt4])
    airflow, fuelflow = convert_mass([airflow, fuelflow])
    pitch_diam *= .0254
    tip_diam *=.0254
    max_diam_casing *= .0254
    min_diam_casing *= .0254
    wall_angle *= np.pi/180
    ref_vel *= .0254*12
    dome_vel_max *= .0254*12
    passage_vel *= .0254*12
    height_turbine_inlet *= .0254

    result = combustor(Tt31, Pt31, airflow, ref_vel, pitch_diam, flow_split, passage_vel, min_diam_casing, max_diam_casing, dome_vel_max, comblendomeheight, fuelflow, LHV, lengthheight, wall_angle, height_turbine_inlet)
    print(result[0])
    print(result[1:3]/.0254)
    print(result[3:5]/.0254)
    print(result[5:7]/.0254)
    print(result[7:]/.0254)

def SonicAreaRatio(M, gamma=1.4):
    '''Returns A*/A from incident Mach number'''
    top = (gamma+1)/2
    bottom = 1+(gamma-1)/2*M**2
    AstarA = M*(top/bottom)**((gamma+1)/(2*(gamma-1)))
    return AstarA

def turbine_vel_diagrams(Tt0, Pt0, Pt2, Pambient, mach2a, Cv, work, alpha2, massflow, tip_diam, hub_diam, spool_speed, stage_eff, gamma_hot, gamma_cold=1.4, R_gas=287.05):
    '''Calculations are acrosss one stage of a turbine. The subscripts denote stations across the stage.
    Station 0 is before the stator, station 1 is between the stator and the rotor, station 2 is after the rotor.
    c = cz + ctj = axial flow + azimuthal flow.
    c is absolute velocity.
    w is realtive to the rotor.'''
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)
    gminus = gamma_hot - 1
    gplus  = gamma_hot + 1
    gminusg = gminus/gamma_hot
    
    pitch_diam1 = (tip_diam + hub_diam)/2
    spool_speed1 = spool_speed*pitch_diam1/2
    # Before stator
    flow_area = np.pi*(tip_diam**2 - hub_diam**2)/4
    densityt0 = Pt0/Tt0/R_gas
    mach0 = rootfind.newton2(find_mach, .3, massflow=massflow, densityt=densityt0, Tt=Tt0, area=flow_area, gamma=gamma_hot, R=R_gas)
    print(mach0)
    Ps0 = isenf.static_pressure(mach0, Pt0, gamma=gamma_hot)
    Ts0 = isenf.static_temperature(mach0, Tt0, gamma=gamma_hot)
    velocity0 = mach0*np.sqrt(gamma_hot*R_gas*Ts0)
    Ps2 = isenf.static_pressure(mach2a, Pt2, gamma=gamma_hot)

    # After stator
    velocity1u = work/spool_speed1
    mach1 = rootfind.newton2(find_mach2, .9, massflow=massflow, densityt=densityt0, Tt=Tt0, flow_area=flow_area, velocity_comp=velocity1u, gamma=gamma_hot, R_gas=R_gas)
    unused, velocity1z = find_mach2(mach1, massflow, densityt0, Tt0, flow_area, velocity1u, gamma=gamma_hot, R_gas=R_gas)
    alpha1 = np.arctan(velocity1u/velocity1z)
    Ps1 = isenf.static_pressure(mach1, Pt0)
    Ts1 = isenf.static_temperature(mach1, Tt0)
    # Before rotor
    reaction = (Ps1**gminusg - Ps2**gminusg)/(Pt0**gminusg - Ps2**gminusg)
    relvel1u = velocity1u - spool_speed1
    beta1 = np.tan(relvel1u/velocity1z)
    relvel1 = np.sqrt(relvel1u**2 + velocity1z**2)
    Ttr1 = Ts1 + relvel1**2/2/cp_hot
    Ptr1 = Ps1*(Ttr1/Ts1)**(1/gminusg)
    # After rotor
    Ttr2 = Ttr1
    machr2 = isenf.temperature_ratio2mach(Ps2/Ptr1, gamma=gamma_hot)
    Ts2i = isenf.static_temperature(machr2, Ttr2, gamma=gamma_hot)
    relvel2i = machr2*np.sqrt(gamma_hot*R_gas*Ts2i)
    relvel2a = np.sqrt(stage_eff)*relvel2i
    Ts2a = Ttr2 - relvel2a**2/2/cp_hot
    densitys2 = Ps2/R_gas/Ts2a
    beta2 = np.arccos(massflow/densitys2/relvel2a/flow_area)
    relvel2u = relvel2a*np.sin(-beta2)
    velocity2z = relvel2a*np.cos(-beta2)
    velocity2u = relvel2u + spool_speed1
    alpha2 = np.tan(velocity2u/velocity2z)
    mach2 = np.sqrt(velocity2u*2+velocity2z**2)/np.sqrt(gamma_hot*R_gas*Ts2a)
    Pt2 = isenf.total_pressure(mach2, Ps2, gamma=gamma_hot)
    Tt2 = isenf.total_temperature(mach2, Ts2a, gamma=gamma_hot)
    print(work, spool_speed1*(velocity1u-velocity2u))
    
    AstarA = SonicAreaRatio(mach2, gamma=gamma_hot)
    Astar = flow_area*AstarA

    PsPt_exit = Pambient/Pt2
    mach_exit = isenf.pressure_ratio2mach(PsPt_exit, gamma_hot)
    velocity_exit = mach_exit*np.sqrt(gamma_hot*R_gas*isenf.total_temperature(mach_exit, Tt2))
    velocity_exit = Cv*velocity_exit
    mach_exit = rootfind.newton2(find_mach3, mach_exit, Tt=Tt2, velocity=velocity_exit, gamma=gamma_hot, R_gas=R_gas)
    AstarAexit = SonicAreaRatio(mach_exit, gamma_hot)
    Aexit = AstarAexit**-1*Astar

    return (reaction, (Tt0, Ts0, Ts1, Tt2, Ts2a), (Pt0, Ps0, Ps1, Pt2, Ps2), (alpha1, beta1, alpha2, beta2), (velocity0, velocity1u, velocity1z, velocity2z), (relvel1u, velocity1z, relvel2u, velocity2z), (mach0, mach1, mach2), (Astar, Aexit))
    
    
def assignment10():
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
    Tt4 = convert_temps(Tt4)
    Pt4 = convert_pressures(Pt4)
    Ps0 = convert_pressures(Ps0)
    Pt49 = convert_pressures(Pt49)
    massflow4 = convert_mass(massflow4)
    work = convert_energy(work)
    tip_diam = tip_diam *.0254
    hub_diam = hub_diam *.0254
    N *= 2*np.pi/60

    res = turbine_vel_diagrams(Tt4, Pt4, Pt49, Ps0, mach49, nozzle_coeff, work, alpha2, massflow4, tip_diam, hub_diam, N, eff_turb, gamma_hot, gamma_cold=1.4, R_gas=287.05)
    print('reaction', np.array([res[0]]))
    print('temp',np.array([res[1]])*9/5)
    print('press',np.array([res[2]])/6895)
    print('angle',np.array([res[3]])*180/np.pi)
    print('abs vel',np.array([res[4]])/.0254/12)
    print('rel vel',np.array([res[5]])/.0254/12)
    print('mach',np.array([res[6]]))
    print('area',np.array([res[7]])/(.0254*12)**2)


if __name__ == '__main__':
    assignment7()
