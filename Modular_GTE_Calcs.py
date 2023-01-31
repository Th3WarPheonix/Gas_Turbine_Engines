
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf

def ambient_properties(mach:float, temperature:float, pressure:float, gamma:float=1.4):
    total_temperature = isenf.T0(mach, temperature, gamma)
    total_pressure    = isenf.p0(mach, pressure, gamma)

    return total_temperature, total_pressure

def inlet(mach:float, total_temperature:float, total_pressure:float, Ipr:float, gamma:float=1.4, R_air:float=287.05):
    total_pressure2 = Ipr*total_pressure
    total_temperature2 = total_temperature
    static_temperature2 = isenf.T(mach, total_temperature2, gamma)
    static_pressure2    = isenf.p(mach, total_pressure2, gamma)
    velocity = mach*np.sqrt(gamma*R_air*static_temperature2)

    return total_temperature2, total_pressure2, static_temperature2, static_pressure2, velocity

def compressor(total_temperature:float, total_pressure:float, efficiency:float, pressure_ratio:float, bypass_ratio:float, gamma:float=1.4, cp_air:float=1006):
    """This set of equations works for any compressor, including a fan. The work results are normalized to the mass flow passing through the compressor."""
    total_pressure2 = pressure_ratio*total_pressure
    total_temperature2_ideal = (total_pressure2/total_pressure)*((gamma-1)/gamma)*total_temperature
    ideal_work = (1+bypass_ratio)*cp_air*(total_temperature2_ideal-total_temperature)
    actual_work = ideal_work/efficiency
    total_temperature2_actual = (total_temperature2_ideal-total_temperature)/efficiency+total_temperature

    return total_temperature2_ideal, total_temperature2_actual, total_pressure2, ideal_work, actual_work
    

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
    Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 1)



if __name__ == '__main__':
    main()