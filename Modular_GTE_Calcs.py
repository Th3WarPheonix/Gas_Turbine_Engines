
import numpy as np
import matplotlib.pyplot as plt
import IsentropicFlow as isenf
import pandas as pd

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

def ambient_properties(mach:float, static_temperature0:float, static_pressure0:float, gamma:float=1.4, R_gas=287.05):
    total_temperature0 = isenf.T0(mach, static_temperature0, gamma)
    total_pressure0    = isenf.p0(mach, static_pressure0, gamma)
    static_density0 = static_pressure0/R_gas/static_temperature0
    total_density0 = total_pressure0/R_gas/total_temperature0
    return total_temperature0, total_pressure0, static_density0, total_density0

def inlet(mach1:float, total_temperature0:float, total_pressure0:float, Ipr:float, gamma:float=1.4, R_gas:float=287.05):
    """Station 1 - Fan Face"""
    total_pressure1     = Ipr*total_pressure0
    total_temperature1  = total_temperature0
    static_temperature1 = isenf.T(mach1, total_temperature1, gamma)
    static_pressure1    = isenf.p(mach1, total_pressure1, gamma)
    static_density1 = static_pressure1/R_gas/static_temperature1
    total_density1  = total_pressure1/R_gas/total_temperature1
    velocity = mach1*np.sqrt(gamma*R_gas*static_temperature1)

    return total_temperature1, total_pressure1, static_temperature1, static_pressure1, velocity, static_density1, total_density1

def compressor(total_temperature:float, total_pressure:float, efficiency:float, pressure_ratio:float, bypass_ratio:float, gamma:float=1.4, R_gas:float=287.05):
    """Station 2, 13 - Fan Exit\n
    Station 3 - Compressor Exit\n
    This set of equations works for any compressor (bypass_ratio=0), including a fan. The work results are normalized to the mass flow passing through the compressor."""
    cp_air = gamma*R_gas/(gamma-1)

    total_pressure2 = pressure_ratio*total_pressure
    total_temperature2_ideal = (total_pressure2/total_pressure)**((gamma-1)/gamma) * total_temperature
    ideal_work = (1+bypass_ratio)*cp_air*(total_temperature2_ideal-total_temperature) # with respect to compressor mass flow
    actual_work = ideal_work/efficiency # with respect to compressor mass flow
    total_temperature2_actual = (total_temperature2_ideal-total_temperature)/efficiency+total_temperature

    return total_temperature2_ideal, total_temperature2_actual, total_pressure2, ideal_work, actual_work

def combustor(Tt31, Tt4, m31, LHV, comb_eff, comb_press_drop, Pt3, gamma_hot, gamma_cold=1.4, R_gas=287.05):
    """Station 4 - Combustor Exit"""
    cp_cold = gamma_cold*R_gas/(gamma_cold-1)
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)
    
    # mfuel = (m31*cp_cold*Tt31 - m31*cp_hot*Tt4)/(cp_hot*Tt4 - LHV) # Does not take into account combustor efficieny, i.e efficiency = 1
    fuel_air_ratio = (cp_cold*Tt4 - cp_cold*Tt31)/(LHV*comb_eff - cp_cold*Tt4)*m31 # Does take into account efficiency mf/m2

    Pt4 = Pt3*comb_press_drop
    return (fuel_air_ratio, Pt4)

def turbine(Tt4, Pt4, Tt31, massflow31, comp_work, fan_work, coolflow, turb_comp_eff, turb_fan_eff, gamma_hot=4/3, gamma=1.4, R_gas=287.05):
    """Station 4.9 - High Pressure Turbine (HPT) Exit\n
    Station 4.95 - Low Pressure Turbine (LPT) Entrance\n
    Station 5 - Low Pressure Turbine Exit"""
    cp_cold = gamma*R_gas/(gamma-1)
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)

    massflow5 = massflow31 + coolflow

    Tt49  = Tt4 - comp_work/(massflow31*cp_hot)
    Pt49 = Pt4*(Tt49/Tt4)**(gamma_hot/(gamma_hot-1))*turb_comp_eff

    Tt495 = (massflow31*cp_hot*Tt49 + coolflow*cp_cold*Tt31)/(massflow5*cp_hot)

    Tt5 = Tt495 - fan_work/(massflow5*cp_hot)
    Pt5 = Pt49*(Tt5/Tt495)**(gamma_hot/(gamma_hot-1))*turb_fan_eff

    return (Tt49, Tt5, Pt49, Pt5, massflow5, Tt495)

def nozzle(Tt5, Pt5, Ps9, Cv, gamma_hot=4/3, gamma=1.4 , R_gas=287.05):
    """Station 9-Nozzle Exit"""
    cp_hot = gamma_hot*R_gas/(gamma_hot-1)

    Ts9i = Tt5*(Ps9/Pt5)**((gamma_hot-1)/gamma_hot)
    vel9i = np.sqrt(2*cp_hot*(Tt5 - Ts9i))
    vel9a = Cv*vel9i
    Ts9a = Tt5 - (vel9a**2)/2/cp_hot
    Pt9 = (Tt5/Ts9i)**(gamma_hot/(gamma_hot-1))*Ps9

    return (Ts9i, Ts9a, Pt9, vel9a)

def engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio, m31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, vel_coeff_core, vel_coeff_fan, thrust, gamma_hot, gamma=1.4, R_gas=287.05):
    if bypass_ratio != 0: # Turbofan
        Tt0, Pt0, rhos0, rhot0   = ambient_properties(mach0, Tambient, Pambient)
        Tt1, Pt1, Ts1, Ps1, Vel1, rhos1, rhot1 = inlet(mach1, Tt0, Pt0, inlet_press_rec)
        Tt13i, Tt13a, Pt13, Wfi, Wfa = compressor(Tt1, Pt1, fan_eff, fan_press_ratio, bypass_ratio) # Fan
        Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt13a, Pt13, comp_eff, comp_press_ratio, 0)

        fuel_air_ratio, Pt4       = combustor(Tt3a, Tt4, m31, LHV, comb_eff, comb_press_drop, Pt3, gamma_hot)
        Tt49, Tt5, Pt49, Pt5, massflow5, Tt495 = turbine(Tt4, Pt4, Tt3a, m31+fuel_air_ratio, Wca, Wfa, turbine_cool_flow, core_turb_eff, fan_turb_eff, gamma_hot)
        Ts9i, Ts9a, Pt9, Vel9     = nozzle(Tt5, Pt5, Pambient, vel_coeff_core) # Core
        Ts19i, Ts19a, Pt19, Vel19 = nozzle(Tt13a, Pt13, Pambient, vel_coeff_fan) # Fan
    else: # Turbojet
        Tt0, Pt0, rhos0, rhot0    = ambient_properties(mach0, Tambient, Pambient)
        Tt1, Pt1, Ts1, Ps1, Vel1, rhos1, rhot1  = inlet(mach1, Tt0, Pt0, inlet_press_rec)
        Tt3i, Tt3a, Pt3, Wci, Wca = compressor(Tt1, Pt1, comp_eff, comp_press_ratio, 0)
        fuel_air_ratio, Pt4       = combustor(Tt3a, Tt4, m31, LHV, comb_eff, comb_press_drop, Pt3, gamma_hot)
        Tt49, Tt5, Pt49, Pt5, massflow5, Tt495 = turbine(Tt4, Pt4, Tt3a, m31+fuel_air_ratio, Wca, 0, turbine_cool_flow, core_turb_eff, fan_turb_eff, gamma_hot)
        Ts9i, Ts9a, Pt9, Vel9     = nozzle(Tt5, Pt5, Pambient, vel_coeff_core) # Core
        Vel19, Tt13i, Tt13a, Wfi, Wfa, Ts19i, Ts19a, Pt13, Pt19 = 0,0,0,0,0,0,0,0,0 # Fan specific variables

    # Calculate the mass flow of air into the compressor based on thrust required
    Vel0 = mach0*np.sqrt(gamma*R_gas*Tambient)
    momentum = massflow5*Vel9 + bypass_ratio*Vel19 - (1+bypass_ratio)*Vel0 # (m/s)/s
    # momentum = massflow5*2762*.0254*12 + bypass_ratio*1184*12*.0254 - (1+bypass_ratio)*Vel0 # (m/s)/s
    massflow2 = thrust*4.4/momentum # kg
    # Calculate the fuel burned base on massflow required
    fuel_density_emp = 6.532 # lbm/gal
    mission_length = 110 # minutes
    fuelflow = massflow2*fuel_air_ratio/m31 # kg/s
    fuelflow_emp = convert_mass(fuelflow, 'imp') # lbm/s
    fuelvol_emp = fuelflow_emp/fuel_density_emp # gal/s
    fuelvol_burned = fuelvol_emp*60*mission_length
    # Calculate the diameter needed to suck in the required mass flow
    totalmassflow = massflow2*(1+bypass_ratio)
    inlet_area = totalmassflow/rhos1/Vel1
    inlet_diameter = np.sqrt(inlet_area*4/np.pi)
    inlet_diameter_emp = inlet_diameter/.0254
    # Conversions
    temperautres_emps = convert_temps([Tambient, Tt0, Tt1, Ts1, Tt13i, Tt13a, Tt3i, Tt3a, Tt4, Tt49, Tt495, Tt5, Ts9i, Ts9a, Ts19i, Ts19a], 'imp')
    Ts0_emp, Tt0_emp, Tt1_emp, Ts1_emp, Tt13i_emp, Tt13a_emp, Tt3i_emp, Tt3a_emp, Tt4_emp, Tt49_emp, Tt495_emp, Tt5_emp, Ts9i_emp, Ts9a_emp, Ts19i_emp, Ts19a_emp = temperautres_emps
    pressures_emps = convert_pressures([Pambient, Pt0, Pt1, Ps1, Pt13, Pt3, Pt4, Pt49, Pt5, Pambient, Pt9, Pt19], 'imp')
    Ps0_emp, Pt0_emp, Pt1_emp, Ps1_emp, Pt13_emp, Pt3_emp, Pt4_emp, Pt49_emp, Pt5_emp, Ps9_emp, Pt9_emp, Pt19_emp = pressures_emps
    Vel1_emp = Vel1 / 12 / .0254
    Vel9_emp = Vel9 / 12 / .0254
    Vel19_emp = Vel19 / 12 / .0254
    Wfi_emp, Wfa_emp, Wci_emp, Wca_emp = convert_work([Wfi, Wfa, Wci, Wca], 'imp')
    massflow2_emp = convert_mass(massflow2, 'imp')
    # Preparing presentation
    labels0 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Mass flow (m/m2)'])
    station0 = pd.Series([Tt0_emp, Ts0_emp, Pt0_emp, Ps0_emp, mach0, 1+bypass_ratio], index=[np.repeat('0 Freestream', len(labels0)), labels0])

    labels1 = np.array(['Total Temperature (R)', 'Static Temperature (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mach', 'Velocity (ft/s)', 'Mass flow (m/m2)'])
    station1 = pd.Series([Tt1_emp, Ts1_emp, Pt1_emp, Ps1_emp, mach1, Vel1_emp, 1+bypass_ratio], index=[np.repeat('1 Fan Inlet', len(labels1)), labels1])

    labels13 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)', 'Mass flow (m/m2)'])
    station13 = pd.Series([Tt13i_emp, Tt13a_emp, Pt13_emp, Wfi_emp, Wfa_emp, bypass_ratio], index=[np.repeat('13 Fan Exit', len(labels13)), labels13])

    labels3 = np.array(['Total Temperature Ideal (R)', 'Total Temperature Actual (R)', 'Total Pressure (psia)', 'Work Ideal (BTU/lbm)', 'Work Actual (BTU/lbm)', 'Mass flow (m/m2)'])
    station3 = pd.Series([Tt3i_emp, Tt3a_emp, Pt3_emp, Wci_emp, Wca_emp, m31], index=[np.repeat('3 Compressor Exit', len(labels3)), labels3])

    labels4 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Fuel/Air', 'Mass flow (m/m2)'])
    station4 = pd.Series([Tt4_emp, Pt4_emp, fuel_air_ratio, m31], index=[np.repeat('4 Combustor Exit', len(labels4)), labels4])

    labels49 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Mass flow (m/m2)'])
    station49 = pd.Series([Tt49_emp, Pt49_emp, m31], index=[np.repeat('4.9 HPT Exit', len(labels49)), labels49])

    labels495 = np.array(['Total Temperature (R)'])
    station495 = pd.Series([Tt495_emp], index=[np.repeat('4.95 LPT Entrance', len(labels495)), labels495])

    labels5 = np.array(['Total Temperature (R)', 'Total Pressure (psia)', 'Mass flow (m/m2)'])
    station5 = pd.Series([Tt5_emp, Pt5_emp, massflow5], index=[np.repeat('5 LPT Exit', len(labels5)), labels5])

    labels9 = np.array(['Static Temperature Ideal (R)', 'Static Temperature Actual (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mass flow (m/m2)', 'Velocity (ft/s)'])
    station9 = pd.Series([Ts9i_emp, Ts9a_emp, Pt9_emp, Ps9_emp, massflow5, Vel9_emp], index=[np.repeat('9 Core Nozzle Exit', len(labels9)), labels9])

    labels19 = np.array(['Total Temperature (R)', 'Static Temperature Ideal (R)', 'Static Temperature Actual (R)', 'Total Pressure (psia)', 'Static Pressure (psia)', 'Mass flow (m/m2)', 'Velocity (ft/s)'])
    station19 = pd.Series([Tt13a_emp, Ts19i_emp, Ts19a_emp, Pt19_emp, Ps9_emp, bypass_ratio, Vel19_emp], index=[np.repeat('19 Fan Nozzle Exit', len(labels19)), labels19])

    labelsSummary = np.array(['Thrust (lbf)', 'Mass Flow 2 (lbm/s)', 'Specific Thrust (lbf/lbm2)', 'Fuel Burn (gal)', 'Inlet Diameter (in)'])
    stationSummary = pd.Series([thrust, massflow2_emp, thrust/massflow2_emp, fuelvol_burned, inlet_diameter_emp], index=[np.repeat('Summary', len(labelsSummary)), labelsSummary])
    
    # Pt0_emp, (1+bypass_ratio)*massflow2_emp
    total_temps = (Tt0_emp, Tt1_emp, Tt13a_emp, Tt3a_emp, Tt3a_emp, Tt4_emp, Tt49_emp, Tt495_emp, Tt5_emp, Tt5_emp, Tt13a_emp, Tt13a_emp)
    total_press = (Pt0_emp, Pt1_emp, Pt13_emp, Pt3_emp, Pt3_emp, Pt4_emp, Pt49_emp, Pt49_emp, Pt5_emp, Pt5_emp, Pt13_emp, Pt13_emp)
    total_mass = ((1+bypass_ratio)*massflow2_emp, (1+bypass_ratio)*massflow2_emp, massflow2_emp, massflow2_emp, massflow2_emp*m31, massflow2_emp*m31, massflow2_emp*m31, massflow2_emp*massflow5, massflow2_emp*massflow5, massflow2_emp*massflow5, bypass_ratio*massflow2_emp, bypass_ratio*massflow2_emp)
    config_temps = pd.Series(total_temps, index=[np.array([0, 1, 2, 3, 3.1, 4, 4.9, 4.95, 5, 9, 13, 19], dtype=str)], name='Temperature (R)')
    config_press = pd.Series(total_press, index=[np.array([0, 1, 2, 3, 3.1, 4, 4.9, 4.95, 5, 9, 13, 19], dtype=str)], name='Pressure (psia)')
    config_mass = pd.Series(total_mass, index=[np.array([0, 1, 2, 3, 3.1, 4, 4.9, 4.95, 5, 9, 13, 19], dtype=str)], name='m (lbm/s)')
    dfResult = pd.DataFrame(pd.concat([config_temps, config_press, config_mass], axis=1))

    return (station0, station1, station3, station4, station49, station495, station5, station9, station13, station19, stationSummary), dfResult

def engine_configurations(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio, comp_eff, comp_press_ratio, m31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, thrust, gamma_hot):
    dfConfigs = pd.DataFrame()
    for i in range(len(comp_press_ratio)):
        stations, dfResult  = engine_walkthrough(Tambient, Pambient, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass_ratio[i], comp_eff, comp_press_ratio[i], m31, LHV, Tt4[i], comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, thrust[i], gamma_hot)
        dfConfigi = pd.DataFrame(pd.concat(stations), columns=['Config {}'.format(i+1)])
        dfConfigs = pd.concat([dfConfigs, dfConfigi], axis=1)
        dfResult.to_csv('Result_Config_{}_{}.csv'.format(comp_press_ratio[i], round(convert_temps(Tt4[i], 'imp')-460, 0)))
    dfConfigs.index = dfConfigs.index.rename(['Station','Property'])
    dfConfigs.to_csv('Engine Configurations.csv')
    return dfConfigs
    
def engine_config_plots(dfConfigs, comp_press_ratio, Tt3max, Tt495max, Fuel_Vol, max_diam):
    fig, axs = plt.subplots(2,2, figsize=(10, 7))
    plt.subplots_adjust(wspace=.4, hspace=.7)

    comp_temp = dfConfigs.loc['3 Compressor Exit', 'Total Temperature Actual (R)'] - 460
    axs[0,0].plot([comp_press_ratio[1], comp_press_ratio[4]], [comp_temp[1], comp_temp[4]], color='green')
    axs[0,0].plot([comp_press_ratio[2], comp_press_ratio[3]], [comp_temp[2], comp_temp[3]], color='green')
    axs[0,0].scatter(comp_press_ratio[1:], comp_temp[1:], marker='.', label='Configs', s=75)
    axs[0,0].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[0,0].set_ylabel('Compressor Exit\nTemperature $T_{t_3}$ ($^\circ$F)')
    axs[0,0].axhline(Tt3max, linestyle='--', label='RFP Limit', color='orange')
    axs[0,0].set_title('Compressor Exit Temperature vs\nCompressor Pressure Ratio')
    axs[0,0].legend()

    Tt495s = convert_temps(dfConfigs.loc['4.95 LPT Entrance', 'Total Temperature (R)'], 'SI')-273
    axs[0,1].plot([comp_press_ratio[1], comp_press_ratio[4]], [Tt495s[1], Tt495s[4]], color='green')
    axs[0,1].plot([comp_press_ratio[2], comp_press_ratio[3]], [Tt495s[2], Tt495s[3]], color='green')
    axs[0,1].scatter(comp_press_ratio[1:], Tt495s[1:], marker='.', label='Configs', s=75)
    axs[0,1].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[0,1].set_ylabel('Exhaust Gas\nTemperature $T_{t_{4.95}}$ ($^\circ$C)')
    axs[0,1].axhline(Tt495max, linestyle='--', label='RFP Limit', color='orange')
    axs[0,1].set_title('Exhaust Gas Temperature vs\nCompressor Pressure Ratio')
    axs[0,1].legend()

    fuel_burn = dfConfigs.loc['Summary', 'Fuel Burn (gal)']
    axs[1,0].plot([comp_press_ratio[1], comp_press_ratio[4]], [fuel_burn[1], fuel_burn[4]], color='green')
    axs[1,0].plot([comp_press_ratio[2], comp_press_ratio[3]], [fuel_burn[2], fuel_burn[3]], color='green')
    axs[1,0].scatter(comp_press_ratio[1:], fuel_burn[1:], marker='.', label='Configs', s=75)
    axs[1,0].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[1,0].set_ylabel('Fuel Burned (gal)')
    axs[1,0].axhline(Fuel_Vol, linestyle='--', label='RFP Limit', color='orange')
    axs[1,0].set_title('Mission fuel burn vs\nCompressor Pressure Ratio')
    axs[1,0].legend()

    diam = dfConfigs.loc['Summary', 'Inlet Diameter (in)']
    axs[1,1].plot([comp_press_ratio[1], comp_press_ratio[4]], [diam[1], diam[4]], color='green')
    axs[1,1].plot([comp_press_ratio[2], comp_press_ratio[3]], [diam[2], diam[3]], color='green')
    axs[1,1].scatter(comp_press_ratio[1:], diam[1:], marker='.', label='Configs', s=75)
    axs[1,1].set_xlabel('Compressor Pressure Ratio $\dfrac{P_{t_3}}{P_{t_2}}$')
    axs[1,1].set_ylabel('Engine Inlet Diameter (in)')
    axs[1,1].axhline(max_diam, linestyle='--', label='RFP Limit', color='orange')
    axs[1,1].set_title('Engine Diameter vs\nCompressor Pressure Ratio')
    axs[1,1].legend()

    plt.savefig('Config_plots.png')
    plt.show()

def rfp2():
    # General values
    alt = 30000 # ft
    thrust = [4800, 4800, 4800, 4800, 4800] # lbf
    spec_trust = 98.4 # lbf/lbmass flow through compressor
    fuel_vol = 1600 # gallons
    # Ambient values
    mach0 = .8
    Ts0 = -44.4 # C
    Ps0 = 4.36 # psia
    # Inlet values
    inlet_press_rec = .99
    # Fan values
    mach1 = .4
    bypass = [2, 0, 0, 0, 0]
    fan_press_ratio = 1.6 # pt13/pt1
    inlet_diam = 51.6 # in
    fan_eff = .88
    # Compressor values
    comp_press_ratio = [7, 7, 9, 7, 9] # pt3/pt2
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
    turbine_inlet_temp = [2100, 2100, 2200, 2200, 2100] # F
    comb_eff = .995
    gamma_hot = 4/3
    # Turbine values
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
    Ts0 += 273
    Ps0 = convert_pressures(Ps0, 'SI')
    Tt4 = convert_temps(np.array(turbine_inlet_temp) + 460, 'SI')
    LHV = convert_work(LHV, 'SI')

    dfConfigs = engine_configurations(Ts0, Ps0, mach0, mach1, inlet_press_rec, fan_eff, fan_press_ratio, bypass, comp_eff, comp_press_ratio, massflow31, LHV, Tt4, comb_eff, comb_press_drop, core_turb_eff, fan_turb_eff, turbine_cool_flow, core_exh_coeff, fan_exh_coeff, thrust, gamma_hot)
    print(dfConfigs)
    engine_config_plots(dfConfigs, comp_press_ratio, max_Tt3, max_Tt495, fuel_vol, max_diam)

    Scaledturbofan = dfConfigs.loc['0 Freestream', 'Total Temperature (R)']

if __name__ == '__main__':
    rfp2()