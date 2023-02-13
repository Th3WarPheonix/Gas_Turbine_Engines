
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf
import pandas as pd

def convert_temps(temps, to):
    """Convert the temperature from K to R"""
    if to == 'emp':
        try:
            temps *= 1.8
            return temps
        except:
            return np.array(temps)*1.8
    elif to == 'SI':
        try:
            temps = temps/1.8
            return temps
        except:
            return np.array(temps)/1.8

def convert_pressures(pressures, to):
    """Convert the pressure from Pa to psia"""
    if to == 'emp':
        try:
            pressures = pressures/6895
            return pressures
        except:
            return np.array(pressures)/6895
    elif to == 'SI':
        try:
            pressures = pressures*6895.0
            return pressures
        except:
            return np.array(pressures)*6895.0

def convert_work(works, to):
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
    
    mfuel = (m31*cp_cold*Tt31 - m31*cp_hot*Tt4)/(cp_hot*Tt4 - LHV)
    fuel_air_ratio = (cp_hot*Tt4 - cp_cold*Tt31)/(LHV*comb_eff - cp_hot*Tt4)
    
    return fuel_air_ratio

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

def compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, Tt31, Pt31, massflow31, mach31, gamma=1.4, R_air=287.05):
    total_area = np.pi*max_tip_diam**2/4 # total cross sectional area of the compressor
    spool_speed = (max_tip_speed)/(max_tip_diam/2) # rad/s
    spool_speed_rpm = spool_speed*(1/(2*np.pi))*(60/1) # rpm
    
    densityt2 = Pt2/(R_air*Tt2)
    Ts31 = isenf.T(mach31, Tt31)
    velocity31 = np.sqrt(gamma*R_air*Ts31)
    densityt31 = Pt31/(R_air*Tt31)
    densitys31 = isenf.density(mach31, densityt31)

    inlet_blade_height = max_tip_diam/(2*(1+inlet_radius_ratio))
    inlet_hub_radius = inlet_blade_height*inlet_radius_ratio
    inlet_hub_area = np.pi*inlet_hub_radius**2
    inlet_flow_area = total_area - inlet_hub_area

    outlet_flow_area = massflow31/densitys31/velocity31
    outlet_hub_area = total_area - outlet_flow_area
    outlet_hub_diam = np.sqrt(outlet_hub_area*4/np.pi)
    outlet_blade_height = (max_tip_diam - outlet_hub_diam)/2

    avg_pitch_diam2  = (max_tip_diam + inlet_hub_radius*2)/2
    avg_pitch_diam31 = (max_tip_diam + outlet_hub_diam)/2
    avg_pitch_diam = (avg_pitch_diam2 + avg_pitch_diam31)/2
    avg_blade_height = (inlet_blade_height + outlet_blade_height)/2
    avg_velocity = spool_speed * avg_pitch_diam/2
    avg_blade_width = avg_blade_height/aspect_ratio
    avg_gap = .25*avg_blade_width
    compressor_length = 6*(2*avg_blade_width+avg_gap) + 5*avg_gap # 6 rotors, 6 stators, 11 gaps in between all of the rotors and stators

    stage_work = work_coeff*avg_velocity**2/2 # work per stage
    num_stages = total_work/stage_work
    num_stages = np.ceil(num_stages)

    mach2 = np.sqrt(2/(gamma-1)*((densityt2*inlet_flow_area*gamma*R_air*Tt2/massflow2)**((gamma-1)/gamma)-1)) # should be between .5 and .55 for conditions in assignment
    Ts2 = isenf.T(mach2, Tt2)
    Ps2 = isenf.p(mach2, Pt2)
    densitys2 = isenf.density(mach2, densityt2)
    velocity2 = mach2*np.sqrt(gamma*R_air*Ts2)

    return ((mach2, Ts2, Ps2, densitys2, velocity2), spool_speed_rpm, num_stages, compressor_length)

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

if __name__ == '__main__':
    Ts3max = 450 # F
    dFrame = engine_configurations()
    engine_config_plots(dFrame, [7, 9, 9], Ts3max, 0, 0)
    # print(assignment5(dFrame))
    # print(assignment6())