
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf
import pandas as pd
import RootFinding as rootfind

def convert_temps(temps, to:str):
    """Convert the temperature from K to R"""
    con_factor = 1.8
    if to == 'imp':
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
    if to == 'imp':
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
    if to == 'imp':
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
    if to == 'imp':
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

def ambient_properties(mach:float, static_temperature0:float, static_pressure0:float, gamma:float=1.4):
    total_temperature0 = isenf.T0(mach, static_temperature0, gamma)
    total_pressure0    = isenf.p0(mach, static_pressure0, gamma)

    return total_temperature0, total_pressure0

def inlet(mach1:float, total_temperature0:float, total_pressure0:float, Ipr:float, gamma:float=1.4, R_gas:float=287.05):
    """Station 1 - Fan Face"""
    total_pressure1     = Ipr*total_pressure0
    total_temperature1  = total_temperature0
    static_temperature1 = isenf.T(mach1, total_temperature1, gamma)
    static_pressure1    = isenf.p(mach1, total_pressure1, gamma)
    velocity = mach1*np.sqrt(gamma*R_gas*static_temperature1)

    return total_temperature1, total_pressure1, static_temperature1, static_pressure1, velocity

def compressor(total_temperature:float, total_pressure:float, efficiency:float, pressure_ratio:float, bypass_ratio:float, gamma:float=1.4, R_gas:float=287.05):
    """Station 2, 13 - Fan Exit\n
    Station 3 - Compressor Exit\n
    This set of equations works for any compressor (bypass_ratio=0), including a fan. The work results are normalized to the mass flow passing through the compressor."""
    cp_air = gamma*R_gas/(gamma-1)
    total_pressure2 = pressure_ratio*total_pressure
    total_temperature2_ideal = (total_pressure2/total_pressure)**((gamma-1)/gamma) * total_temperature
    ideal_work = (1+bypass_ratio) * cp_air * (total_temperature2_ideal-total_temperature) # work with respect to compressor mass flow
    actual_work = ideal_work/efficiency # work with respect to compressor mass flow
    total_temperature2_actual = (total_temperature2_ideal-total_temperature)/efficiency+total_temperature

    return total_temperature2_ideal, total_temperature2_actual, total_pressure2, ideal_work, actual_work

def combustor(Tt31, Tt4, m31, LHV, comb_eff, comb_press_drop, Pt3, gamma_hot, gamma_cold=1.4, R_gas=287.05):
    """Station 4 - Combustor Exit"""
    cp_cold = gamma_cold*R_gas/(gamma_cold-1)
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)
    
    mfuel = (m31*cp_cold*Tt31 - m31*cp_hot*Tt4)/(cp_hot*Tt4 - LHV) # Does not take into account combustor efficieny, i.e efficiency = 1
    fuel_air_ratio = (cp_hot*Tt4 - cp_cold*Tt31)/(LHV*comb_eff - cp_hot*Tt4) # Does take into account efficiency
    
    Pt4 = Pt3*comb_press_drop
    return (fuel_air_ratio, Pt4)

def turbine(Tt4, Pt4, massflow31, comp_work, fan_work, coolflow, turb_comp_eff, turb_fan_eff, gamma_hot=1.3, gamma=1.4, R_gas=287.05):
    """Station 4.9 - High Pressure Turbine (HPT) Exit\n
    Station 4.95 - Low Pressure Turbine (LPT) Entrance\n
    Station 5 - Low Pressure Turbine Exit"""
    cp_cold = gamma*R_gas/(gamma-1)
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)

    massflow5 = massflow31 + coolflow

    Tt49  = (massflow31*cp_hot*Tt4 - comp_work/turb_comp_eff)/(massflow31*cp_hot)
    Tt5 = (massflow5*cp_hot*Tt49 - fan_work/turb_fan_eff)/(massflow5*cp_hot)
    Pt49 = Pt4*(Tt49/Tt4)**(gamma_hot/(gamma_hot-1))
    Pt5 = Pt49*(Tt5/Tt49)**(gamma_hot/(gamma_hot-1))

    return (Tt49, Tt5, Pt49, Pt5, massflow5)

def nozzle(Tt5, Pt5, Ps9, densitys9, Cv, gamma_hot=1.3, gamma=1.4 , R_gas=287.05):
    """Station 9-Nozzle Exit"""
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)
    Ts9i = (Tt5/(Pt5/Ps9))**((gamma_hot-1)/gamma_hot)
    vel9i = 2*np.sqrt(cp_hot*(Tt5 - Ts9i))
    vel9a = Cv*vel9i
    Ts9a = cp_hot*Tt5 - vel9a**2/2
    Pt9 = Ps9 + 1/2*densitys9*vel9a**2
    return (Ts9i, Ts9a, Pt9, vel9a)

def inlet_design(stream_density:float, stream_velocity:float, massflow:float, A0AHL:float, mach_throat:float, Tt0:float, Pt0:float, Ts1:float, Ps1:float, mach_fan:float, diffuser_angle:float, gamma:float=1.4, R_gas:float=287.05):

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

    '''Fan face calcs'''
    fan_static_temp  = isenf.T(mach_fan, Tt0)
    fan_static_press = isenf.p(mach_fan, Pt0)
    speed_of_sound = np.sqrt(gamma*R_gas*fan_static_temp)
    fan_velocity = speed_of_sound*mach_fan
    fan_density = fan_static_press/(R_gas*fan_static_temp)
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
    return (massflow/densitys/mach/sound_speed - area,)

def find_c2z(c2z, w2t, Ptr2, Ttr2, massflow, flow_area2, gamma=1.4, R_gas=287.05):
    cp_air = gamma*R_gas/(gamma-1)
    w2z = c2z
    w2 = w2z + w2t*1j
    Ts2 = Ttr2 - np.abs(w2)**2/2/cp_air
    Ps2 = Ptr2*(Ts2/Ttr2)**(gamma/(gamma-1))
    densitys2 = Ps2/Ts2/R_gas
    cz2 = massflow/densitys2/flow_area2
    return (cz2 - cz2,)

def find_efficiency(stage_eff, c1, Tt1, Pt1, Pt3_req, spool_speed1, spool_speed2, Ptr2, Ttr2, massflow, flow_area2, loss_coeffs, gamma=1.4, R_gas=287.05):
    cp_air = gamma*R_gas/(gamma-1)
    delta_ct = ((Pt3_req/Pt1)**((gamma-1)/gamma)-1)*cp_air*Tt1/spool_speed1/stage_eff # change in c theta across the rotor
        
    c2z = np.real(c1)
    c2t = np.imag(c1) + delta_ct
    w2t = c2t - spool_speed2
    # ans = rootfind.newton2(find_c2z, x1=133, w2t=w2t, Ptr2=Ptr2, Ttr2=Ttr2, massflow=massflow, flow_area2=flow_area2)
    ans = 133
    c2 = ans + c2t*1j
    w2 = c2 - spool_speed2*1j

    Ts2 = Ttr2 - np.abs(w2)**2/2/cp_air
    Ps2 = Ptr2*(Ts2/Ttr2)**(gamma/(gamma-1))
    densitys2 = Ps2/Ts2/R_gas
    Tt2 = Ts2 + np.abs(c2)**2/2/cp_air
    Pt2 = Ps2*(Ts2/Tt2)**-(gamma/(gamma-1))  
    Pt3 = Pt2 - loss_coeffs*densitys2*np.abs(c2)**2/2
    return Pt3_req - Pt3, c2, w2

def compressor_design(max_tip_diam, max_tip_speed, aspect_ratio, work_coeff, total_work, inlet_radius_ratio, Tt2, Pt2, massflow2, Tt31, Pt31, massflow31, mach31, gamma=1.4, R_gas=287.05):
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

    mach2 = rootfind.newton2(find_area, .469, massflow=massflow2, densityt=densityt2, Tt=Tt2, area=inlet_flow_area)
    Ts2 = isenf.T(mach2, Tt2)
    Ps2 = isenf.p(mach2, Pt2)
    densitys2 = isenf.density(mach2, densityt2)
    velocity2 = mach2*np.sqrt(gamma*R_gas*Ts2)
    
    return ((inlet_hub_radius*2, outlet_hub_diam, avg_gap, avg_blade_height), spool_speed_rpm, num_stages, compressor_length)

def engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio, m31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, vel_coeff_core, vel_coeff_fan, gamma_hot, R_gas=287.05):
    Tt0, Pt0 = ambient_properties(mach0, Tambient, Pambient)
    Tt1, Pt1, Ts1, Ps1, Vel1 = inlet(mach1, Tt0, Pt0, inlet_press_rec)
    Tt13i, Tt13a, Pt13, Wfi, Wfa = compressor(Tt1, Pt1, fan_eff, fan_press_ratio, bypass_ratio) # for the fan
    Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 0)
    Tt31 = Tt13a
    Wc = Wca
    Wf = Wfa
    fuel_air_ratio, Pt4 = combustor(Tt31, Tt4, m31, LHV, comb_eff, comb_press_drop, Pt3, gamma_hot)
    Tt49, Tt5, Pt49, Pt5, massflow5 = turbine(Tt4, Pt4, m31, Wc, Wf, turbine_cool_flow, core_turb_eff, fan_turb_eff, gamma_hot)
    Ts9i, Ts9a, Pt9, Vel9 = nozzle(Tt5, Pt5, Pambient, Pambient/R_gas/Tambient, vel_coeff_core)

    Thrust = 7000
    momentum = massflow5*Vel9 - bypass_ratio*Vel1
    massflow2 = Thrust*4.4/momentum

    temperautres_emps = convert_temps([Tambient, Tt0, Tt1, Ts1, Tt13i, Tt13a, Tt3i, Tt3a, Tt4, Tt49, Tt5, Ts9i, Ts9a], 'imp')
    Ts0_emp, Tt0_emp, Tt1_emp, Ts1_emp, Tt13i_emp, Tt13a_emp, Tt3i_emp, Tt3a_emp, Tt4_emp, Tt49_emp, Tt5_emp, Ts9i_emp, Ts9a_emp = temperautres_emps
    pressures_emps = convert_pressures([Pambient, Pt0, Pt1, Ps1, Pt13, Pt3, Pt4, Pt49, Pt5, Pambient, Pt9], 'imp')
    Ps0_emp, Pt0_emp, Pt1_emp, Ps1_emp, Pt13_emp, Pt3_emp, Pt4_emp, Pt49_emp, Pt5_emp, Ps9_emp, Pt9_emp = pressures_emps
    Vel1_emp = Vel1 / 12 / .0254
    Vel9_emp = Vel9 / 12 / .0254
    Wfi_emp, Wfa_emp, Wci_emp, Wca_emp = convert_work([Wfi, Wfa, Wci, Wca], 'imp')
    massflow2_emp = convert_mass(massflow2, 'imp')

    labels0 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Mass flow (m/m2)'])
    station0 = pd.Series([Tt0_emp, Ts0_emp, Pt0_emp, Ps0_emp, mach0, bypass_ratio], index=[np.repeat('0 Freestream', len(labels0)), labels0])

    labels1 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Velocity (ft/s)', 'Mass flow (m/m2)'])
    station1 = pd.Series([Tt1_emp, Ts1_emp, Pt1_emp, Ps1_emp, mach1, Vel1_emp, bypass_ratio], index=[np.repeat('1 Fan Inlet', len(labels1)), labels1])

    labels13 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)', 'Mass flow (m/m2)'])
    station13 = pd.Series([Tt13i_emp, Tt13a_emp, Pt13_emp, Wfi_emp, Wfa_emp, bypass_ratio], index=[np.repeat('13 Bypass/Compressor Inlet', len(labels13)), labels13])

    labels3 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)', 'Mass flow (m/m2)'])
    station3 = pd.Series([Tt3i_emp, Tt3a_emp, Pt3_emp, Wci_emp, Wca_emp, m31], index=[np.repeat('3 Compressor Exit', len(labels3)), labels3])

    labels4 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Fuel/Air', 'Mass flow (m/m2)'])
    station4 = pd.Series([Tt4_emp, Pt4_emp, fuel_air_ratio, m31], index=[np.repeat('4 Combustor Exit', len(labels4)), labels4])

    labels49 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Mass flow (m/m2)'])
    station49 = pd.Series([Tt49_emp, Pt49_emp, m31], index=[np.repeat('4.9 HPT Exit', len(labels49)), labels49])

    labels5 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Mass flow (m/m2)'])
    station5 = pd.Series([Tt5_emp, Pt5_emp, massflow5], index=[np.repeat('5 LPT Exit', len(labels5)), labels5])

    labels9 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mass flow (m/m2)'])
    station9 = pd.Series([Ts9i_emp, Ts9a_emp, Pt9_emp, Ps9_emp, massflow5], index=[np.repeat('9 LPT Exit', len(labels9)), labels9])

    labelsSummary = np.array(['Thrust (lbf)', 'Mass Flow 2 (lbm/s)', 'Specific Thrust (lbf/lbm2)'])
    stationSummary = pd.Series([Thrust, massflow2_emp, Thrust/massflow2], index=[np.repeat('Summary LPT Exit', len(labelsSummary)), labelsSummary])
    
    return (station0, station1, station13, station3, station4, station49, station5, station9, stationSummary)

def engine_configurations(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio, m31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, gamma_hot):
    dfConfigs = pd.DataFrame()
    for i in range(len(comp_press_ratio)):
        stations  = engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio[i], m31, LHV, Tt4[i], comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, gamma_hot)
        dfConfigi = pd.DataFrame(pd.concat(stations), columns=['Config {}'.format(i+1)])
        dfConfigs = pd.concat([dfConfigs, dfConfigi], axis=1)
    dfConfigs.index = dfConfigs.index.rename(['Station','Property'])
    dfConfigs.to_csv('Engine Configurations.csv')
    return dfConfigs
    
def engine_config_plots(dfConfigs, comp_press_ratio, Tt3max, Tt495max):
    fig, axs = plt.subplots(2,2, figsize=(10, 6))
    plt.subplots_adjust(wspace=.4, hspace=.4)

    axs[0,0].plot(comp_press_ratio, dfConfigs.loc['3 Compressor Exit', 'Total Temperature Actual (R)'] - 460)
    axs[0,0].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[0,0].set_ylabel('Compressor Exit\nTemperature $T_{t_3}$ ($^\circ$F)')
    axs[0,0].axhline(Tt3max, linestyle=':')
    axs[0,0].set_title('Compressor Exit Temperature vs\nCompressor Pressure Ratio')

    axs[0,1].plot(comp_press_ratio, dfConfigs.loc['4.9 HPT Exit', 'Total Temperature (R)'] - 460)
    axs[0,1].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[0,1].set_ylabel('Exhaust Gas\nTemperature $T_{t_{4.95}}$ ($^\circ$F)')
    axs[0,1].axhline(Tt495max, linestyle=':')
    axs[0,1].set_title('Exhaust Gas Temperature vs\nCompressor Pressure Ratio')

    axs[2,0].plot(comp_press_ratio, dfConfigs.loc['4.9 HPT Exit', 'Total Temperature (R)'] - 460)
    axs[2,0].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[2,0].set_ylabel('Exhaust Gas\nTemperature $T_{t_{4.95}}$ ($^\circ$F)')
    axs[2,0].axhline(Tt495max, linestyle=':')
    axs[2,0].set_title('Mission fuel burn vs\nCompressor Pressure Ratio')

    # axs[2,1].plot(comp_press_ratio, dfConfigs.loc['4.9 HPT Exit', 'Total Temperature (R)'] - 460)
    # axs[2,1].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    # axs[2,1].set_ylabel('Exhaust Gas\nTemperature $T_{t_{4.95}}$ ($^\circ$F)')
    # axs[2,1].axhline(Tt495max, linestyle=':')
    # axs[2,1].set_title('Engine Diameter vs\nCompressor Pressure Ratio')

    plt.show()

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

def compressor_vel_diagrams(Tt1, Pt1, massflow, alpha1, press_ratio, num_stages, Dt1, Dp1, Dp2, flow_area3, spool_speed, stage_eff, loss_coeffr, loss_coeffs, reaction, alpha3, gamma=1.4, R_gas=287.05):
    '''Calculations are acrosss one stage of a compressor. The subscripts denote stations across the stage.
    Station 1 is before the rotor, station 2 is between the rotor and the stator, station 3 is after the stator.
    c = cz + ctj = axial flow + azimuthal flow.
    c is absolute velocity.
    w is realtive to the rotor.'''
    cp_air = gamma*R_gas/(gamma-1)
    spool_speed1 = spool_speed*Dp1/2
    spool_speed2 = spool_speed*Dp2/2
    
    hub_diam1 = 2*Dp1 - Dt1
    hub_area1 = np.pi*hub_diam1**2/4
    flow_area1 = np.pi*Dt1**2/4 - hub_area1

    hub_diam2 = 2*Dp2 - Dt1
    hub_area2 = np.pi*hub_diam2**2/4
    flow_area2 = np.pi*Dt1**2/4 - hub_area2

    densityt1 = Pt1/R_gas/Tt1
    mach1 = rootfind.newton2(find_area, .5, massflow=massflow, densityt=densityt1, Tt=Tt1, area=flow_area1)
    Ts1 = isenf.T(mach1, Tt1)
    Ps1 = isenf.p(mach1, Pt1)
    densitys1 = isenf.density(mach1, densityt1)
    sound_speed1 = np.sqrt(gamma*R_gas*Ts1)
    air_vel1 = mach1*sound_speed1
    
    c1 = air_vel1*np.cos(alpha1) + air_vel1*np.sin(alpha1)*1j
    w1 = c1 - spool_speed1*1j
    beta1 = np.arctan(np.imag(w1)/np.real(w1))

    Ttr1 = Ts1 + np.abs(w1)**2/2/cp_air
    Ptr1 = Ps1*(Ttr1/Ts1)**(gamma/(gamma-1))

    Ttr2 = Ttr1 - spool_speed1**2/2/cp_air + spool_speed2**2/2/cp_air
    Ptr2 = Ptr1 - loss_coeffr*densitys1*np.abs(w1)**2/2
   
    Pt3_req = Pt1 * press_ratio**(1/num_stages)
    stage_eff = rootfind.newton2(find_efficiency, stage_eff, c1=c1, Tt1=Tt1, Pt1=Pt1, Pt3_req=Pt3_req, spool_speed1=spool_speed1, spool_speed2=spool_speed2, Ptr2=Ptr2, Ttr2=Ttr2, massflow=massflow, flow_area2=flow_area2, loss_coeffs=loss_coeffs)
    unused, c2, w2 = find_efficiency(stage_eff, c1, Tt1, Pt1, Pt3_req, spool_speed1, spool_speed2, Ptr2, Ttr2, massflow, flow_area2, loss_coeffs)
    alpha2 = np.arctan(np.imag(c2)/np.real(c2))
    beta2 = np.arctan(np.imag(w2)/np.real(w2))
    delta_ct = ((Pt3_req/Pt1)**((gamma-1)/gamma)-1)*cp_air*Tt1/spool_speed1/stage_eff
    work_stage = spool_speed1*delta_ct # specific work per stage

    Ts2 = Ttr2 - np.abs(w2)**2/2/cp_air
    Ps2 = Ptr2*(Ts2/Ttr2)**(gamma/(gamma-1))
    densitys2 = Ps2/Ts2/R_gas
    Tt2 = Ts2 + np.abs(c2)**2/2/cp_air
    Pt2 = Ps2*(Ts2/Tt2)**-(gamma/(gamma-1))
    Pt3 = Pt2 - loss_coeffs*densitys2*np.abs(c2)**2/2

    Tt3 = Tt2
    reaction = (np.abs(w1)**2-np.abs(w2)**2)/2/spool_speed2/delta_ct
    densityt3 = Pt3/R_gas/Tt3
    mach3 = rootfind.newton2(find_area, .5, massflow=massflow, densityt=densityt3, Tt=Tt3, area=flow_area3)
    Ts3 = isenf.T(mach3, Tt3)
    Ps3 = isenf.p(mach3, Pt3)
    densitys3 = isenf.density(mach3, densityt3)
    sound_speed3 = np.sqrt(gamma*R_gas*Ts3)
    air_vel3 = mach3*sound_speed3
    c3 = np.cos(alpha3)*air_vel3 + np.sin(alpha3)*air_vel3 
          
    print()
    print('beta1 deg', beta1*180/np.pi)
    print('C1 ft/s', c1/.0254/12)
    print('W1 ft/s', w1/.0254/12)
    print('Ptr1 psia', convert_pressures(Ptr1, 'imp'))
    print('Ttr1 R', convert_temps(Ttr1, 'imp'))
    print('Wheel speed ft/s', spool_speed1/.0254/12)
    print('||w1|| ft/s', abs(w1)/.0254/12)
    print()
    print('beta2 deg', beta2*180/np.pi)
    print('C2 ft/s', c2/.0254/12)
    print('W2 ft/s', w2/.0254/12)
    print('Ptr2 psia', convert_pressures(Ptr2, 'imp'))
    print('Ttr2 R', convert_temps(Ttr2, 'imp'))
    print('||w2|| ft/s', abs(w2)/.0254/12)
    print('||c2|| ft/s', abs(c2)/.0254/12)
    print('alpha2 deg', alpha2*180/np.pi)
    print('Pt2 psia', convert_pressures(Pt2, 'imp'))
    print('Tt2 R', convert_temps(Tt2, 'imp'))
    print('Ps2 psia', convert_pressures(Ps2, 'imp'))
    print('Ts2 R', convert_temps(Ts2, 'imp'))
    print('Stage work btu/lbm', convert_work(work_stage, 'imp'))
    print()
    print('Pt3 psia', convert_pressures(Pt3, 'imp'))
    print('Tt3 R', convert_temps(Tt3, 'imp'))
    print('reaction', reaction)
    print('stage efficiency', stage_eff)
    print('c3', c3/.0254/12)
    print('Ps3 psia', convert_pressures(Ps3, 'imp'))

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
        
    compressor_vel_diagrams(Tt1, Pt1, massflow1, alpha1, comp_press_ratio, num_stages, Dt1, Dp1, Dp2, area3, spool_speed_rads, stage_eff, loss_coeff_r, loss_coeff_s, reaction, alpha3)

def rfp2():
    # General values
    alt = 30000 # ft
    thrust = 7000 # lbf
    spec_trust = 98.4 # lbf/lbmass flow through compressor
    # Ambient values
    mach0 = .8
    Ts0 = -44.4 # C
    Ps0 = 4.36 # psia
    # Inlet values
    inlet_press_rec = .99
    # Fan values
    mach1 = .4
    bypass = 2
    fan_press_ratio = 1.6 # pt13/pt1
    inlet_diam = 51.6 # in
    fan_eff = .88
    # Compressor values
    comp_press_ratio = [7]#[7, 9, 9] # pt3/pt2
    massflow2 = 71.2 # lbm
    comp_eff = .87
    comp_leak_flow = .01
    turbine_cool_flow = .06
    customer_bleed_flow = .03
    massflow31 = 1 - comp_leak_flow - turbine_cool_flow - customer_bleed_flow
    # Combustor values
    comb_press_drop = .95 # Pt4/Pt3
    comb_eff = .995 # actual heat/ ideal heat
    LHV = 18550 # BTU/lbmfuel
    turbine_inlet_temp = [2100] #[2100, 2100, 2200] # F
    comb_eff = .995
    gamma_hot = 1.3
    # Turbine values
    Tt4 = 2100 # F
    core_turb_eff = .92
    fan_turb_eff = .925
    # Nozzle values
    core_exh_coeff = .983 # V9actual/V9ideal
    fan_exh_coeff = .985
    # Operating limits
    max_diam = 30 # in
    max_Tt3 = 450 # F
    max_Tt495 = 890 # C
    max_Tt9 = 880 # C
    # Conversions
    Ts0 = Ts0 + 273
    Ps0 = convert_pressures(Ps0, 'SI')
    Tt4 = convert_temps(turbine_inlet_temp, 'SI')
    LHV = convert_work(LHV, 'SI')

    dfConfigs = engine_configurations(Ts0, Ps0, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass, comp_eff, comp_press_ratio, massflow31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, gamma_hot)
    print(dfConfigs)
    # engine_config_plots(dfConfigs, comp_press_ratio, max_Tt3, max_Tt495)

if __name__ == '__main__':
    Ts3max = 450 # F
    # dFrame = engine_configurations()
    # engine_config_plots(dFrame, [7, 9, 9], Ts3max, 0, 0)
    # print(assignment5(dFrame))
    # ans6 = assignment6()
    # assignment7()
    # airfoil_count()
    rfp2()