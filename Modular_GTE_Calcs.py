
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf
import pandas as pd
import RootFinding as rootfind

def convert_temps(temps, to:str):
    """Convert the temperature from K to R"""
    con_factor = 1.8
    if to == 'emp':
        try:
            temps *= con_factor
            return temps
        except:
            return np.array(temps)*con_factor
    elif to == 'SI':
        try:
            temps = temps/con_factor
            return temps
        except:
            return np.array(temps)/con_factor

def convert_mass(mass, to:str):
    """Convert the temperature from lbm to kg"""
    con_factor = 2.20462
    if to == 'emp':
        try:
            mass *= con_factor
            return mass
        except:
            return np.array(mass)*con_factor
    elif to == 'SI':
        try:
            mass = mass/con_factor
            return mass
        except:
            return np.array(mass)/con_factor


def convert_pressures(pressures, to:str):
    """Convert the pressure from Pa to psia"""
    con_factor = 6895.0
    if to == 'emp':
        try:
            pressures = pressures/con_factor
            return pressures
        except:
            return np.array(pressures)/con_factor
    elif to == 'SI':
        try:
            pressures = pressures*con_factor
            return pressures
        except:
            return np.array(pressures)*con_factor

def convert_work(works, to:str):
    """Convert the work from J to Btu"""
    if to == 'emp':
        try:
            works = works / 1055 / 2.205
            return works
        except:
            return np.array(works) / 1055 / 2.205
    elif to == 'SI':
        try:
            works = works * 1055 * 2.205
            return works
        except:
            return np.array(works) * 1055 * 2.205

def ambient_properties(mach:float, temperature:float, pressure:float, gamma:float=1.4):
    total_temperature = isenf.T0(mach, temperature, gamma)
    total_pressure    = isenf.p0(mach, pressure, gamma)

    return total_temperature, total_pressure

def inlet(mach:float, total_temperature:float, total_pressure:float, Ipr:float, gamma:float=1.4, R_air:float=287.05):
    total_pressure2     = Ipr*total_pressure
    total_temperature2  = total_temperature
    static_temperature2 = isenf.T(mach, total_temperature2, gamma)
    static_pressure2    = isenf.p(mach, total_pressure2, gamma)
    velocity = mach*np.sqrt(gamma*R_air*static_temperature2)

    return total_temperature2, total_pressure2, static_temperature2, static_pressure2, velocity

def compressor(total_temperature:float, total_pressure:float, efficiency:float, pressure_ratio:float, bypass_ratio:float, gamma:float=1.4, R_air:float=287.05):
    """This set of equations works for any compressor, including a fan. The work results are normalized to the mass flow passing through the compressor."""
    cp_air = gamma*R_air/(gamma-1)
    total_pressure2 = pressure_ratio*total_pressure
    total_temperature2_ideal = (total_pressure2/total_pressure)**((gamma-1)/gamma) * total_temperature
    ideal_work = (1+bypass_ratio) * cp_air * (total_temperature2_ideal-total_temperature)
    actual_work = ideal_work/efficiency
    total_temperature2_actual = (total_temperature2_ideal-total_temperature)/efficiency+total_temperature

    return total_temperature2_ideal, total_temperature2_actual, total_pressure2, ideal_work, actual_work

def combustor(Tt31, Tt4, m31, LHV, comb_eff, gamma_hot, gamma_cold=1.4, R_air=287.05):
    cp_cold = gamma_cold*R_air/(gamma_cold-1)
    cp_hot = gamma_hot*R_air/(gamma_hot-1)
    
    mfuel = (m31*cp_cold*Tt31 - m31*cp_hot*Tt4)/(cp_hot*Tt4 - LHV) # Does not take into account combustor efficieny, i.e efficiency = 1
    fuel_air_ratio = (cp_hot*Tt4 - cp_cold*Tt31)/(LHV*comb_eff - cp_hot*Tt4) # Does take into account efficiency
    
    return fuel_air_ratio

def nozzle(Tt, Pt, pressure_ratio):
    # efficiency = 
    pass

def inlet_design(stream_density:float, stream_velocity:float, massflow:float, A0AHL:float, mach_throat:float, Tt0:float, Pt0:float, Ts1:float, Ps1:float, mach_fan:float, diffuser_angle:float, gamma:float=1.4, R_air:float=287.05):

    '''Entrance calcs'''
    streamtube_area = massflow/(stream_velocity*stream_density)
    streamtube_diameter = np.sqrt(4*streamtube_area/np.pi)
    highlight_area = streamtube_area/A0AHL
    highlight_diameter = np.sqrt(4*highlight_area/np.pi)

    '''Throat calcs'''
    throat_static_temp  = isenf.T(mach_throat, Tt0)
    throat_static_press = isenf.p(mach_throat, Pt0)
    speed_of_sound = np.sqrt(gamma*R_air*throat_static_temp)
    throat_velocity = speed_of_sound*mach_throat
    throat_density = throat_static_press/(R_air*throat_static_temp)
    throat_area = massflow/(throat_velocity*throat_density)
    throat_diameter = np.sqrt(4*throat_area/np.pi)

    '''Fan face calcs'''
    fan_static_temp  = isenf.T(mach_fan, Tt0)
    fan_static_press = isenf.p(mach_fan, Pt0)
    speed_of_sound = np.sqrt(gamma*R_air*fan_static_temp)
    fan_velocity = speed_of_sound*mach_fan
    fan_density = fan_static_press/(R_air*fan_static_temp)
    fan_area = massflow/(fan_velocity*fan_density)
    fan_diameter = np.sqrt(4*fan_area/np.pi)

    '''Diffuser length calcs'''
    throat_radius = throat_diameter/2
    fan_radius = fan_diameter/2
    diffuser_length = (fan_radius-throat_radius)/np.tan(diffuser_angle*np.pi/180)

    return np.array([streamtube_diameter, highlight_diameter, throat_diameter, fan_diameter]), diffuser_length, diffuser_length/throat_diameter

def find_area(mach, massflow, densityt, Tt, area, gamma=1.4, R=287.05):
    densitys = isenf.density(mach, densityt)
    Ts = isenf.T(mach, Tt)
    sound_speed = np.sqrt(gamma*R*Ts)
    return massflow/densitys/mach/sound_speed - area

def compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, Tt31, Pt31, massflow31, mach31, gamma=1.4, R_air=287.05):
    total_area = np.pi*max_tip_diam**2/4 # total cross sectional area of the compressor
    spool_speed = (max_tip_speed)/(max_tip_diam/2) # rad/s
    spool_speed_rpm = spool_speed*(1/(2*np.pi))*(60/1) # rpm
    gap_width_ratio = .25 # gap width to blade width ratio

    densityt2 = Pt2/(R_air*Tt2)

    Ts31 = isenf.T(mach31, Tt31)

    velocity31 = mach31*np.sqrt(gamma*R_air*Ts31)
    densityt31 = Pt31/(R_air*Tt31)
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

    mach2 = rootfind.newton2(find_area, .469, massflow=massflow2, densityt=densityt2, Tt=Tt2, area=inlet_flow_area)
    Ts2 = isenf.T(mach2, Tt2)
    Ps2 = isenf.p(mach2, Pt2)
    densitys2 = isenf.density(mach2, densityt2)
    velocity2 = mach2*np.sqrt(gamma*R_air*Ts2)
    
    return ((inlet_hub_radius*2, outlet_hub_diam, avg_gap, avg_blade_height), spool_speed_rpm, num_stages, compressor_length)

def engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio, m31, LHV, Tt4, comb_eff, gamma_hot):
    Tt0, Pt0 = ambient_properties(mach0, Tambient, Pambient)
    Tt1, Pt1, Ts1, Ps1, Vel1 = inlet(mach1, Tt0, Pt0, inlet_press_rec)
    Tt13i, Tt13a, Pt13, Wfi, Wfa = compressor(Tt1, Pt1, fan_eff, fan_press_ratio, bypass_ratio) # for the fan
    Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 0)
    Tt31 = Tt13a
    fuel_air_ratio = combustor(Tt31, Tt4, m31, LHV, comb_eff, gamma_hot)

    temperautres_emps = convert_temps([Tambient, Tt0, Tt1, Ts1, Tt13i, Tt13a, Tt3i, Tt3a], 'emp')
    Ts0_emp, Tt0_emp, Tt1_emp, Ts1_emp, Tt13i_emp, Tt13a_emp, Tt3i_emp, Tt3a_emp = temperautres_emps

    pressures_emps = convert_pressures([Pambient, Pt0, Pt1, Ps1, Pt13, Pt3], 'emp')
    Ps0_emp, Pt0_emp, Pt1_emp, Ps1_emp, Pt13_emp, Pt3_emp = pressures_emps

    Vel1_emp = Vel1 / 12 / .0254

    works_emps = convert_work([Wfi, Wfa, Wci, Wca], 'emp')
    Wfi_emp, Wfa_emp, Wci_emp, Wca_emp = works_emps

    labels0 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach'])
    station0 = pd.Series([Tt0_emp, Ts0_emp, Pt0_emp, Ps0_emp, mach0], index=[np.repeat('0 Freestream', len(labels0)), labels0])

    labels1 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Velocity (ft/s)'])
    station1 = pd.Series([Tt1_emp, Ts1_emp, Pt1_emp, Ps1_emp, mach1, Vel1_emp], index=[np.repeat('1 Fan Inlet', len(labels1)), labels1])

    labels13 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station13 = pd.Series([Tt13i_emp, Tt13a_emp, Pt13_emp, Wfi_emp, Wfa_emp], index=[np.repeat('13 Bypass/Compressor Inlet', len(labels13)), labels13])

    labels3 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station3 = pd.Series([Tt3i_emp, Tt3a_emp, Pt3_emp, Wci_emp, Wca_emp], index=[np.repeat('3 Compressor Exit', len(labels3)), labels3])
    
    return (station0, station1, station13, station3)

def engine_configurations():
    Rair = 287.05 # J/kg-K
    mach0 = .8
    Tambient = -44.4 # C
    Pambient = 4.36 # psia

    Tambient = Tambient + 273 # K
    Pambient = Pambient * 6894.76 # Pa

    bypass_ratio = 2
    inlet_press_rec = .99
    mach1 = .4

    fan_eff = .88
    fan_press_ratio = 1.6
    
    comp_eff = .87
    comp_press_ratio = [7, 9, 9]

    m31 = .9
    LHV = 18550 # BTU/lbmfuel
    LHV = convert_work(LHV, 'SI')
    turbine_inlet_temp = [2100, 2100, 2200] # F
    Tt4 = convert_temps(turbine_inlet_temp, 'SI')
    comb_eff = .995
    gamma_hot = 1.3

    dfConfigs = pd.DataFrame()
    for i in range(3):
        station0, station1, station13, station3 = engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio[i], m31, LHV, Tt4[i], comb_eff, gamma_hot)
        dfConfigi = pd.DataFrame(pd.concat([station0, station1, station13, station3]), columns=['Config {}'.format(i+1)])
        dfConfigs = pd.concat([dfConfigs, dfConfigi], axis=1)
    dfConfigs.index = dfConfigs.index.rename(['Station','Property'])
    dfConfigs.to_csv('Engine Configurations.csv')
    return dfConfigs
    
def engine_config_plots(dfConfigs, comp_press_ratio, Ts3max, T495max, Tt9max):
    fig, axs = plt.subplots(2,2)
    axs[0,0].plot(comp_press_ratio, dfConfigs.loc['3 Compressor Exit', 'Total Temperature Actual (R)'] - 460)
    axs[0,0].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[0,0].set_ylabel('Compressor Exit\nTemperature $T_{s_3}$ ($^\circ$F)')
    axs[0,0].axhline(Ts3max, linestyle=':')
    axs[0,0].set_title('Compressor Pressure Ratio vs\nCompressor Exit Temperature')
    plt.show()

def assignment5(dfConfigs, R_air=287.05):
    Tt0 = convert_temps(dfConfigs['Config 1'].loc['0 Freestream', 'Total Temperature (R)'], 'SI')
    Pt0 = convert_pressures(dfConfigs['Config 1'].loc['0 Freestream', 'Total Pressure (psia)'], 'SI')
    Ts0 = convert_temps(dfConfigs['Config 1'].loc['0 Freestream', 'Static Temperature (R)'], 'SI')
    Ps0 = convert_pressures(dfConfigs['Config 1'].loc['0 Freestream', 'Static Pressure (psia)'], 'SI')
    Ts1 = convert_temps(dfConfigs['Config 1'].loc['1 Fan Inlet', 'Static Temperature (R)'], 'SI')
    Ps1 = convert_pressures(dfConfigs['Config 1'].loc['1 Fan Inlet', 'Static Pressure (psia)'], 'SI')

    '''Inlet Design'''
    density0 = Ps0/(R_air*Ts0)
    velocity0 = 795.5 # ft/s
    velocity0 = velocity0 * 12*.0254
    massflow0 = 70.17 # lbm/s
    massflow0 = massflow0/2.205
    A0Ahl = .7
    throat_mach = .7
    fan_mach = .4
    diffuser_angle = 5
    diams, diff_length, diff_LH = inlet_design(density0, velocity0, massflow0, A0Ahl, throat_mach, Tt0, Pt0, Ts1, Ps1, fan_mach, diffuser_angle)
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
    total_work = convert_work(total_work, 'SI') # J/(kg/s)

    return compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, Tt31, Pt31, massflow31, mach31)

def compressor_vel_diagrams(Tt1, Pt1, massflow, alpha1, press_ratio, num_stages, Dt1, Dp1, Dp2, flow_area3, spool_speed, stage_eff, loss_coeffr, loss_coeffs, reaction, gamma=1.4, R_air=287.05):
    '''Calculations are acrosss one stage of a compressor. The subscripts denote stations across the stage.
    Station 1 is before the rotor, station 2 is between the rotor and the stator, station 3 is after the stator.
    c = cz + ctj = axial flow + azimuthal flow.
    c is absolute velocity.
    w is realtive to the rotor.'''
    cp_air = gamma*R_air/(gamma-1)
    spool_speed1 = spool_speed*Dp1/2
    spool_speed2 = spool_speed*Dp2/2
    
    hub_diam1 = 2*Dp1 - Dt1
    hub_area1 = np.pi*hub_diam1**2/4
    flow_area1 = np.pi*Dt1**2/4 - hub_area1

    hub_diam2 = 2*Dp2 - Dt1
    hub_area2 = np.pi*hub_diam2**2/4
    flow_area2 = np.pi*Dt1**2/4 - hub_area2

    densityt1 = Pt1/R_air/Tt1
    mach1 = rootfind.newton2(find_area, .5, massflow=massflow, densityt=densityt1, Tt=Tt1, area=flow_area1)
    Ts1 = isenf.T(mach1, Tt1)
    Ps1 = isenf.p(mach1, Pt1)
    densitys1 = isenf.density(mach1, densityt1)
    sound_speed1 = np.sqrt(gamma*R_air*Ts1)
    air_vel1 = mach1*sound_speed1
    
    c1 = air_vel1*np.cos(alpha1) + air_vel1*np.sin(alpha1)*1j
    w1 = c1 - spool_speed1*1j
    beta1 = np.arctan(np.imag(w1)/np.real(w1))

    Ttr1 = Ts1 + np.abs(w1)**2/2/cp_air
    Ttr2 = Ttr1 - spool_speed1**2/2/cp_air + spool_speed2**2/2/cp_air

    beta2 = np.arctan(-2*(reaction-0.5)*spool_speed1/np.real(c1)-np.tan(alpha1))
    w2t = -2*reaction*spool_speed1 - np.imag(w1)
    w2z = w2t/np.tan(beta2)
    w2 = w2z + w2z*1j
    c2 = w2 + spool_speed1*1j
    alpha2 = np.arctan(np.imag(c2)/ np.real(c2))
    c2_mag = np.abs(c2)
    Pt2 = Pt1 - .5*loss_coeffr*densitys1*np.abs(w1)**2
    
    machz = np.real(c1)/sound_speed1
    macht = spool_speed1/sound_speed1
    bottom = 1/macht**2 + (gamma-1)/(2*np.cos(alpha1)*np.cos(alpha1))*(machz/macht)**2
    first = (gamma-1)/bottom
    second = 1 + (machz/macht)*(np.tan(beta2)-np.tan(alpha1))
    Tt2 = (1 + first*second)*Tt1

    Pt3 = Pt1 * press_ratio**(1/num_stages)
    # Pt3 = Pt2 - .5*loss_coeffs*densitys2*np.abs(c2)**2
    Tt3 = ((Pt3/Pt1)**((gamma-1)/gamma)-1)*Tt1 + Tt1
    densityt3 = Pt3/R_air/Tt3
    mach3 = rootfind.newton2(find_area, .5, massflow=massflow, densityt=densityt3, Tt=Tt3, area=flow_area3)
    Ts3 = isenf.T(mach3, Tt3)
    Ps3 = isenf.p(mach3, Pt3)
    sound_speed3 = np.sqrt(gamma*R_air*Ts3)
    air_vel3 = mach3*sound_speed3

    delta_ct = ((Pt3/Pt1)**((gamma-1)/gamma)-1)*cp_air*Tt1/spool_speed1 # change in c theta across the rotor
    work_stage = spool_speed1*delta_ct*massflow

    Ts2 = reaction*(Tt3-Tt1) + Ts1

    fig, (axs1, axs2) = plt.subplots(1, 2)
    axs1.axhline(y=0)
    axs1.plot([0,np.real(c1)], [0, np.imag(c1)], label='c1', color='blue')
    axs1.plot([0,np.real(w1)], [0, np.imag(w1)], label='w1', color='red')
    axs1.plot([np.real(c1), np.real(c1)], [np.imag(c1), np.imag(c1)-spool_speed1], label='Spool1', color='green')
    axs1.legend()

    axs2.axhline(0)
    axs2.plot([0,np.real(c2)], [0, np.imag(c2)], label='c2', color='blue')
    axs2.plot([0,np.real(w2)], [0, np.imag(w2)], label='w2', color='red')
    axs2.plot([np.real(w2), np.real(w2)], [np.imag(w2), np.imag(w2)+spool_speed1], label='Spool2', color='green')
    axs2.legend()
    plt.show()

def airfoil_count():
    rotor_solidity = 1.3
    rotor_pitch_diam = 20.75 * .0254
    rotor_airfoil_width = 3.56 *.0254

    stator_solidity = 1.2
    stator_pitch_diam = 22*.0254
    stator_airfoil_width = 3.3 *.0254

    meanline_slope = 5.9 # deg

    beta1 = 59.4*np.pi/180
    beta2 = 43.3*np.pi/180
    alpha2 = 25.2*np.pi/180
    alpha3 = 0*np.pi/180

    rotor_stagger_angle = (beta1 + beta2)/2

    rotor_chord = rotor_airfoil_width/np.cos(rotor_stagger_angle)

    rotor_num_airfoils  = np.pi*rotor_pitch_diam/pitch_spacing
    stator_num_airfoils = np.pi*stator_pitch_diam/pitch_spacing

def assignment7():
    Tt1 = 464.5 # R
    Pt1 = 6.58 # psia
    massflow1 = 70.17 # lbm/s
    Dp1 = 20.71 # in
    Dt1 = 29.58 # in
    alpha1 = 1e-10 # deg

    Dp2 = 21.37 # in

    area3 = 522.1 # in^2
    alpha3 = 0 # deg

    spool_speed_rpm = 11619 # rpm
    reaction = .5
    stage_eff = .87
    comp_press_ratio = 8
    num_stages = 6
    loss_coeff_r = .1
    loss_coeff_s = .07

    Tt1 = convert_temps(Tt1, "SI")
    Pt1 = convert_pressures(Pt1, "SI")
    massflow1 = convert_mass(massflow1, "SI")
    Dp1 = Dp1 * .0254
    Dt1 = Dt1 * .0254
    alpha1 = alpha1 * np.pi/180
    Dp2 = Dp2 * .0254
    area3 = area3*.0254**2
    alpha3 = alpha3 * np.pi/180

    spool_speed_rads = spool_speed_rpm*2*np.pi/60
        
    compressor_vel_diagrams(Tt1, Pt1, massflow1, alpha1, comp_press_ratio, num_stages, Dt1, Dp1, Dp2, area3, spool_speed_rads, stage_eff, loss_coeff_r, loss_coeff_s, reaction)


if __name__ == '__main__':
    Ts3max = 450 # F
    # dFrame = engine_configurations()
    # engine_config_plots(dFrame, [7, 9, 9], Ts3max, 0, 0)
    # print(assignment5(dFrame))
    # ans6 = assignment6()
    assignment7()