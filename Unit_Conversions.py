
import numpy as np

def convert_temperature(temps, units:str='K'):
    """Convert temperature between Kelvin and Rankine
    units = 'K' for converting to Kelvin
    units = 'R' for converting to Rankine"""
    factor = 1.8
    match units:
        case 'R':
            try:
                temps *= factor
                return temps
            except:
                return np.array(temps)*factor
        case 'K':
            try:
                temps = temps/factor
                return temps
            except:
                return np.array(temps)/factor
        case _:
            print('Did not convert temperature')

def convert_mass(mass, units:str='kg'):
    """Convert mass between lbm and kg
    units = 'lbm' for converting to lbm
    units = 'kg' for converting to kg"""
    factor = 2.20462
    match units:
        case 'lbm':
            try:
                mass *= factor
                return mass
            except:
                return np.array(mass)*factor
        case 'kg':
            try:
                mass = mass/factor
                return mass
            except:
                return np.array(mass)/factor
        case _:
            print('Did not convert mass')

def convert_pressure(pressures, units:str='Pa'):
    """Convert pressure between Pascals and PSI
    units = 'Pa' for converting to Pascals
    units = 'psi' for converting to PSI"""
    factor = 6895.0
    match units:
        case 'psi':
            try:
                pressures = pressures/factor
                return pressures
            except:
                return np.array(pressures)/factor
        case 'Pa':
            try:
                pressures = pressures*factor
                return pressures
            except:
                return np.array(pressures)*factor
        case _:
            print('Did not convert pressure')

def convert_energy(energy, units:str='J'):
    """Convert mass specific energy/work between Btu/lbm and J/kg
    units = 'J' for converting to J/kg
    units = 'BTU' for converting to BTU/lbm"""
    factor = 1055 * 2.205
    match units:
        case 'BTU':
            try:
                energy = energy/factor
                return energy
            except:
                return np.array(energy)/factor
        case 'J':
            try:
                energy = energy*factor
                return energy
            except:
                return np.array(energy)*factor
        case _:
            print('Did not convert energy')

def convert_force(force, units:str='N'):
    """Convert force between Newtons and lbf
    units = 'N' for converting to Newtons
    tunitso = 'lbf' for converting to lbf"""
    factor = 4.44822
    match units:
        case 'lbf':
            try:
                force = force/factor
                return force
            except:
                return np.array(force)/factor
        case 'N':
            try:
                force = force*factor
                return force
            except:
                return np.array(force)*factor
        case _:
            print('Did not convert force')

def convert_length(length, units:str='m in'):
    """Convert length from meters to inches or meters to feet NOT between feet and inches
    units = 'm in' for converting to meters from inches
    units = 'm ft' for converting to meters from feet
    units = 'in' for converting to inches
    units = 'ft' for converting to feet"""
    factor1 = .0254
    factor2 = .0254*12
    match units:
        case 'm in':
            try:
                length = length/factor1
                return length
            except:
                return np.array(length)/factor1
        case 'm ft':
            try:
                length = length/factor2
                return length
            except:
                return np.array(length)/factor2
        case 'in':
            try:
                length = length*factor1
                return length
            except:
                return np.array(length)*factor1
        case 'ft':
            try:
                length = length*factor2
                return length
            except:
                return np.array(length)*factor2
        case _:
            print('Did not convert length')

def convert_area(area, units:str='m in'):
    """Convert length from sq meters to sq inches or sq meters to sq feet NOT between sq feet and sq inches
    units = 'm in' for converting to sq meters from sq inches
    units = 'm ft' for converting to sq meters from sq feet
    units = 'in' for converting to sq inches
    units = 'ft' for converting to sq feet"""
    factor1 = .0254**2
    factor2 = (.0254*12)**2
    match units:
        case 'm in':
            try:
                area = area/factor1
                return area
            except:
                return np.array(area)/factor1
        case 'm ft':
            try:
                area = area/factor2
                return area
            except:
                return np.array(area)/factor2
        case 'in':
            try:
                area = area*factor1
                return area
            except:
                return np.array(area)*factor1
        case 'ft':
            try:
                area = area*factor2
                return area
            except:
                return np.array(area)*factor2
        case _:
            print('Did not convert area')

def convert_speed(speed, units:str='m in'):
    """Convert length from meters/sec to inches/sec or meters/sec to feet/sec NOT between feet and inches
    \nunits = 'm in' for converting to meters/sec from inches/sec
    \nunits = 'm ft' for converting to meters/sec from ft/sec
    \nunits = 'in' for converting to inches/sec
    \nunits = 'ft' for converting to feet/sec
    \nunits = 'mph m' for converting to mph to meters/sec
    \nunits = 'kts mph' for converting to knots to mph
    """
    factor1 = .0254**2 # m/s to in/s
    factor2 = .0254*12 # m/s to ft/s
    factor3 = 5280/3600 # mph to ft/s
    factor4 = 5280/3600*12 # mph to in/s
    factor5 = 1.150779 # knots to mph
    match units:
        case 'm in':
            try:
                speed = speed/factor1
                return speed
            except:
                return np.array(speed)/factor1
        case 'm ft':
            try:
                speed = speed/factor2
                return speed
            except:
                return np.array(speed)/factor2
        case 'in':
            try:
                speed = speed*factor1
                return speed
            except:
                return np.array(speed)*factor1
        case 'ft':
            try:
                speed = speed*factor2
                return speed
            except:
                return np.array(speed)*factor2
        case 'mph ft':
            try:
                speed = speed*factor3
                return speed
            except:
                return np.array(speed)*factor3
        case 'mph m':
            try:
                speed = speed*factor3*factor2
                return speed
            except:
                return np.array(speed)*factor3*factor2
        case 'mph in':
            try:
                speed = speed*factor4
                return speed
            except:
                return np.array(speed)*factor4
        case 'ft mph':
            try:
                speed = speed/factor3
                return speed
            except:
                return np.array(speed)/factor3
        case 'm mph':
            try:
                speed = speed/factor3/factor2
                return speed
            except:
                return np.array(speed)/factor3/factor2
        case 'in mph':
            try:
                speed = speed/factor4
                return speed
            except:
                return np.array(speed)/factor4
        case 'mph kts':
            try:
                speed = speed/factor5
                return speed
            except:
                return np.array(speed)/factor5
        case 'kts mph':
            try:
                speed = speed*factor5
                return speed
            except:
                return np.array(speed)*factor5
        case _:
            print('Did not convert speed')


