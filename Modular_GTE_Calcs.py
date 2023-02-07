
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf
import pandas as pd

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

def convert_temps(temps):
    """Convert the temperature from K to R"""
    try:
        temps *= 1.8
        return np.array([temps])
    except:
        return np.array(temps)*1.8

def convert_pressures(pressures):
    """Convert the pressure from Pa to psia"""
    try:
        pressures *= 1.8
        return np.array([pressures])
    except:
        return np.array(pressures)/6895

def convert_work(works):
    """Convert the work from J to Btu"""
    try:
        works *= 1.8
        return np.array([works])
    except:
        return np.array(works) / 1055 / 2.205

def main():
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
    comp_press_ratio = 7
    
    Tt0, Pt0 = ambient_properties(mach0, Tambient, Pambient)
    Tt1, Pt1, Ts1, Ps1, Vel1 = inlet(mach1, Tt0, Pt0, inlet_press_rec)
    Tt13i, Tt13a, Pt13, Wfi, Wfa = compressor(Tt1, Pt1, fan_eff, fan_press_ratio, bypass_ratio) # for the fan
    Tt2i, Tt2a, Pt2, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 0)

    temperautres_emps = convert_temps([Tambient, Tt0, Tt1, Ts1, Tt13i, Tt13a, Tt2i, Tt2a])
    Ts0_emp, Tt0_emp, Tt1_emp, Ts1_emp, Tt13i_emp, Tt13a_emp, Tt3i_emp, Tt3a_emp = temperautres_emps

    pressures_emps = convert_pressures([Pambient, Pt0, Pt1, Ps1, Pt13, Pt2])
    Ps0_emp, Pt0_emp, Pt1_emp, Ps1_emp, Pt13_emp, Pt2_emp = pressures_emps

    Vel1_emp = Vel1 / 12 / .0254

    works_emps = convert_work([Wfi, Wfa, Wci, Wca])
    Wfi_emp, Wfa_emp, Wci_emp, Wca_emp = works_emps

    labels0 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach'])
    station0 = pd.Series([Tt0_emp, Ts0_emp, Pt0_emp, Ps0_emp, mach0], index=[np.repeat('0 Freestream', len(labels0)), labels0])

    labels1 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Velocity (ft/s)'])
    station1 = pd.Series([Tt1_emp, Ts1_emp, Pt1_emp, Ps1_emp, mach1, Vel1_emp], index=[np.repeat('1 Fan Inlet', len(labels1)), labels1])

    labels13 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station13 = pd.Series([Tt13i_emp, Tt13a_emp, Pt13_emp, Wfi_emp, Wfa_emp], index=[np.repeat('13 Bypass Inlet', len(labels13)), labels13])

    labels2 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station2 = pd.Series([Tt3i_emp, Tt3a_emp, Pt2_emp, Wci_emp, Wca_emp], index=[np.repeat('2 Compressor Inlet', len(labels2)), labels2])
    
    dfConfigs = pd.DataFrame(pd.concat([station0, station1, station13, station2]), columns=['Config1'])
    dfConfigs.index=dfConfigs.index.rename(['Station','Property'])
    print(dfConfigs)
    density0 = Pambient/(Rair*Tambient)
    velocity0 = 795.5 # ft/s
    velocity0 = velocity0 * 12*.0254
    massflow0 = 70.17 # lbm/s
    massflow0 = massflow0/2.205
    A0Ahl = .7
    throat_mach = .7
    fan_mach = .4
    diffuser_angle = 5
    diams, diff_length, diff_LH = inlet_design(density0, velocity0, massflow0, A0Ahl, throat_mach, Tt0, Pt0, Ts1, Ps1, fan_mach, diffuser_angle)
    print(diams/.0254, diff_length/.0254, diff_LH)

if __name__ == '__main__':
    main()