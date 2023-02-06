
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

def inlet_design():
    pass

def main():
    mach0 = .8
    Tambient = -44.4 # C
    Pambient = 4.36 # psia

    Tambient = Tambient + 273
    Pambient = Pambient * 6894.76

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
    Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 0)

    Ts0_emp = Tambient * 1.8
    Tt0_emp = Tt0 * 1.8
    Tt1_emp = Tt1 * 1.8
    Ts1_emp = Ts1 * 1.8
    Tt13i_emp = Tt13i * 1.8
    Tt13a_emp = Tt13a * 1.8
    Tt3i_emp = Tt3i * 1.8
    Tt3a_emp = Tt3a * 1.8

    Ps0_emp = Pambient / 6895
    Pt0_emp = Pt0 / 6895
    Pt1_emp = Pt1 / 6895
    Ps1_emp = Ps1 / 6895
    Pt13_emp = Pt13 / 6895
    Pt3_emp = Pt3 / 6895

    Vel1_emp = Vel1 / 12 / .0254

    Wfi_emp = Wfi / 1055 / 2.205
    Wfa_emp = Wfa / 1055 / 2.205
    Wci_emp = Wci / 1055 / 2.205
    Wca_emp = Wca / 1055 / 2.205

    labels0 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach'])
    station0 = pd.Series([Tt0_emp, Ts0_emp, Pt0_emp, Ps0_emp, mach0], index=[np.zeros_like(labels0, dtype=int), labels0])

    labels1 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Velocity (ft/s)'])
    station1 = pd.Series([Tt1_emp, Ts1_emp, Pt1_emp, Ps1_emp, mach1, Vel1_emp], index=[np.ones_like(labels1, dtype=int), labels1])

    labels13 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station13 = pd.Series([Tt13i_emp, Tt13a_emp, Pt13_emp, Wfi_emp, Wfa_emp], index=[np.ones_like(labels13, dtype=int)*13, labels13])

    labels3 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)'])
    station3 = pd.Series([Tt3i_emp, Tt3a_emp, Pt3_emp, Wci_emp, Wca_emp], index=[np.ones_like(labels3, dtype=int)*3, labels3])
    
    dfConfigs = pd.DataFrame(pd.concat([station0, station1, station13, station3]), columns=['Config1'])
    dfConfigs.index=dfConfigs.index.rename(['Station','Property'])
    print(dfConfigs)
    
if __name__ == '__main__':
    main()