from ambiance import Atmosphere
import numpy as np
import matplotlib.pyplot as plt

"""
All inputs and outputs are in SI units. Functions output absolute temperatures"""

class stdday(Atmosphere):

    @property
    def _setup(self):
        rad0 = 6356.577 # km, earth radius
        geopot_alt = rad0*self.h/1000/(rad0+self.h/1000)

        z = [0, 11, 20, 32, 47, 51, 71, 84.852]
        L = [-6.5, 0, 1, 2.8, 0, -2.8, -2, 0]

        layer0 = ((z[0] <= geopot_alt)*(geopot_alt < z[1])).astype(int)
        layer1 = ((z[1] <= geopot_alt)*(geopot_alt < z[2])).astype(int)
        layer2 = ((z[2] <= geopot_alt)*(geopot_alt < z[3])).astype(int)
        layer3 = ((z[3] <= geopot_alt)*(geopot_alt < z[4])).astype(int)
        layer4 = ((z[4] <= geopot_alt)*(geopot_alt < z[5])).astype(int)
        layer5 = ((z[5] <= geopot_alt)*(geopot_alt < z[6])).astype(int)
        layer6 = ((z[6] <= geopot_alt)*(geopot_alt < z[7])).astype(int)
        layer7 = ((z[7] <= geopot_alt)).astype(int)

        T0 = 288.15
        T1 = T0 + L[0]*(z[1]-z[0])
        T2 = T1 + L[1]*(z[2]-z[1])
        T3 = T2 + L[2]*(z[3]-z[2])
        T4 = T3 + L[3]*(z[4]-z[3])
        T5 = T4 + L[4]*(z[5]-z[4])
        T6 = T5 + L[5]*(z[6]-z[5])
        T7 = T6 + L[6]*(z[7]-z[6])

        T = [T0, T1, T2, T3, T4, T5, T6, T7]
        layer=[layer0, layer1, layer2, layer3, layer4, layer5, layer6, layer7]

        return T, layer, geopot_alt, L, z
    
    @property
    def temperature(self):

        T, layer, gpa, L, z = self._setup
        temp = 0    + (T[0] + L[0]*(gpa-z[0]))*layer[0]
        temp = temp + (T[1] + L[1]*(gpa-z[1]))*layer[1]
        temp = temp + (T[2] + L[2]*(gpa-z[2]))*layer[2]
        temp = temp + (T[3] + L[3]*(gpa-z[3]))*layer[3]
        temp = temp + (T[4] + L[4]*(gpa-z[4]))*layer[4]
        temp = temp + (T[5] + L[5]*(gpa-z[5]))*layer[5]
        temp = temp + (T[6] + L[6]*(gpa-z[6]))*layer[6]
        temp = temp + (T[7] + L[7]*(gpa-z[7]))*layer[7]
        return temp
    
    @property
    def pressure(self):
        g0 = 9806.65 # mm/s^2
        Rgas = 287.05 # J/mol-K
               
        T, layer, gpa, L, z = self._setup

        P0 = 1
        P1 = P0*(T[1]/T[0])**(-g0/Rgas/L[0])
        P2 = P1*np.exp(-g0*(z[2]-z[1])/Rgas/T[1])
        P3 = P2*(T[3]/T[2])**(-g0/Rgas/L[2])
        P4 = P3*(T[4]/T[3])**(-g0/Rgas/L[3])
        P5 = P4*np.exp(-g0*(z[5]-z[4])/Rgas/T[4])
        P6 = P5*(T[6]/T[5])**(-g0/Rgas/L[5])

        prss = (self.temperature/T[0])**(-g0/Rgas/L[0])*layer[0]
        prss = prss + P1*np.exp(-g0*(gpa-z[1])/Rgas/T[1])*layer[1]
        prss = prss + P2*(self.temperature/T[2])**(-g0/Rgas/L[2])*layer[2]
        prss = prss + P3*(self.temperature/T[3])**(-g0/Rgas/L[3])*layer[3]
        prss = prss + P4*np.exp(-g0*(gpa-z[4])/Rgas/T[4])*layer[4]
        prss = prss + P5*(self.temperature/T[5])**(-g0/Rgas/L[5])*layer[5]
        prss = prss + P6*(self.temperature/T[6])**(-g0/Rgas/L[6])*layer[6]
        return prss*101325
    

class coldday(Atmosphere):

    @property
    def _setup(self):
        rad0 = 6356.577 # km, earth radius
        geopot_alt = rad0*self.h/1000/(rad0+self.h/1000)

        z = [0, 1, 3, 9.5, 13, 15.5, 18.5, 22.5]
        L = [25, 0, -6, 0, -8.88, 0, 4.6, -.775]

        layer0 = ((z[0] <= geopot_alt)*(geopot_alt < z[1])).astype(int)
        layer1 = ((z[1] <= geopot_alt)*(geopot_alt < z[2])).astype(int)
        layer2 = ((z[2] <= geopot_alt)*(geopot_alt < z[3])).astype(int)
        layer3 = ((z[3] <= geopot_alt)*(geopot_alt < z[4])).astype(int)
        layer4 = ((z[4] <= geopot_alt)*(geopot_alt < z[5])).astype(int)
        layer5 = ((z[5] <= geopot_alt)*(geopot_alt < z[6])).astype(int)
        layer6 = ((z[6] <= geopot_alt)*(geopot_alt < z[7])).astype(int)
        layer7 = ((z[7] <= geopot_alt)).astype(int)
        
        T0 = 222.10
        T1 = T0 + L[0]*(z[1]-z[0])
        T2 = T1 + L[1]*(z[2]-z[1])
        T3 = T2 + L[2]*(z[3]-z[2])
        T4 = T3 + L[3]*(z[4]-z[3])
        T5 = T4 + L[4]*(z[5]-z[4])
        T6 = T5 + L[5]*(z[6]-z[5])
        T7 = T6 + L[6]*(z[7]-z[6])

        T = [T0, T1, T2, T3, T4, T5, T6, T7]
        layer=[layer0, layer1, layer2, layer3, layer4, layer5, layer6, layer7]

        return (T, layer, geopot_alt, L, z)
    
    @property
    def temperature(self):

        T, layer, gpa, L, z = self._setup
        temp = (T[0] + L[0]*(gpa-z[0]))*layer[0]
        temp = temp + (T[1] + L[1]*(gpa-z[1]))*layer[1]
        temp = temp + (T[2] + L[2]*(gpa-z[2]))*layer[2]
        temp = temp + (T[3] + L[3]*(gpa-z[3]))*layer[3]
        temp = temp + (T[4] + L[4]*(gpa-z[4]))*layer[4]
        temp = temp + (T[5] + L[5]*(gpa-z[5]))*layer[5]
        temp = temp + (T[6] + L[6]*(gpa-z[6]))*layer[6]
        temp = temp + (T[7] + L[7]*(gpa-z[7]))*layer[7]
  
        return temp
    
    @property
    def pressure(self):
        g0 = 9806.65 # mm/s^2
        Rgas = 287.05 # J/mol-K
               
        T, layer, gpa, L, z = self._setup

        P0 = 1
        P1 = P0*(T[1]/T[0])**(-g0/Rgas/L[0])
        P2 = P1*np.exp(-g0*(z[2]-z[1])/Rgas/T[1])
        P3 = P2*(T[3]/T[2])**(-g0/Rgas/L[2])
        P4 = P3*np.exp(-g0*(z[4]-z[3])/Rgas/T[3])
        P5 = P4*(T[5]/T[4])**(-g0/Rgas/L[4])
        P6 = P5*np.exp(-g0*(z[6]-z[5])/Rgas/T[5])
        P7 = P6*np.exp(-g0*(z[7]-z[6])/Rgas/T[6])

        prss = (self.temperature/T[0])**(-g0/Rgas/L[0])*layer[0]
        prss = prss + P1*np.exp(-g0*(gpa-z[1])/Rgas/T[1])*layer[1]
        prss = prss + P2*(self.temperature/T[2])**(-g0/Rgas/L[2])*layer[2]
        prss = prss + P3*np.exp(-g0*(gpa-z[3])/Rgas/T[3])*layer[3]
        prss = prss + P4*(self.temperature/T[4])**(-g0/Rgas/L[4])*layer[4]
        prss = prss + P5*np.exp(-g0*(gpa-z[5])/Rgas/T[5])*layer[5]
        prss = prss + P6*(self.temperature/T[6])**(-g0/Rgas/L[6])*layer[6]
        prss = prss + P7*(self.temperature/T[7])**(-g0/Rgas/L[7])*layer[7]
        print(P1, P2, P3, P4, P5, P6, P7)
        return prss*101325

class hotday(Atmosphere):

    @property
    def _setup(self):
        rad0 = 6356.577 # km, earth radius
        geopot_alt = rad0*self.h/1000/(rad0+self.h/1000)

        z = [0, 12, 20.5]
        L = [-7, .8, 1.4]

        layer0 = ((z[0] <= geopot_alt)*(geopot_alt < z[1])).astype(int)
        layer1 = ((z[1] <= geopot_alt)*(geopot_alt < z[2])).astype(int)
        layer2 = ((z[2] <= geopot_alt)).astype(int)
        
        T0 = 312.6
        T1 = T0 + L[0]*(z[1]-z[0])
        T2 = T1 + L[1]*(z[2]-z[1])

        T = [T0, T1, T2]
        layer = [layer0, layer1, layer2]

        return (T, layer, geopot_alt, L, z)
    
    @property
    def temperature(self):

        T, layer, gpa, L, z = self._setup
        temp = (T[0] + L[0]*(gpa-z[0]))*layer[0]
        temp = temp + (T[1] + L[1]*(gpa-z[1]))*layer[1]
        temp = temp + (T[2] + L[2]*(gpa-z[2]))*layer[2]

        return temp
    
    @property
    def pressure(self):
        g0 = 9806.65 # mm/s^2
        Rgas = 287.05 # J/mol-K
               
        T, layer, gpa, L, z = self._setup

        P0 = 1
        P1 = P0*(T[1]/T[0])**(-g0/Rgas/L[0])
        P2 = P1*(T[2]/T[1])**(-g0/Rgas/L[1])

        prss = (self.temperature/T[0])**(-g0/Rgas/L[0])*layer[0]
        prss = prss + P1*np.exp(-g0*(gpa-z[1])/Rgas/T[1])*layer[1]
        prss = prss + P2*np.exp(-g0*(gpa-z[2])/Rgas/T[2])*layer[2]

        return prss*101325

class tropday(Atmosphere):

    @property
    def _setup(self):
        rad0 = 6356.577 # km, earth radius
        geopot_alt = rad0*self.h/1000/(rad0+self.h/1000)

        z = [0, 16, 21]
        L = [-7, 3.8, 2.48]

        layer0 = ((z[0] <= geopot_alt)*(geopot_alt < z[1])).astype(int)
        layer1 = ((z[1] <= geopot_alt)*(geopot_alt < z[2])).astype(int)
        layer2 = ((z[2] <= geopot_alt)).astype(int)
        
        T0 = 305.27
        T1 = T0 + L[0]*(z[1]-z[0])
        T2 = T1 + L[1]*(z[2]-z[1])

        T = [T0, T1, T2]
        layer = [layer0, layer1, layer2]

        return T, layer, geopot_alt, L, z
    
    @property
    def temperature(self):

        T, layer, gpa, L, z = self._setup
        temp = (T[0] + L[0]*(gpa-z[0]))*layer[0]
        temp = temp + (T[1] + L[1]*(gpa-z[1]))*layer[1]
        temp = temp + (T[2] + L[2]*(gpa-z[2]))*layer[2]

        return temp
    
    @property
    def pressure(self):
        g0 = 9806.65 # mm/s^2
        Rgas = 287.05 # J/mol-K
               
        T, layer, gpa, L, z = self._setup

        P0 = 1
        P1 = P0*(T[1]/T[0])**(-g0/Rgas/L[0])
        P2 = P1*(T[2]/T[1])**(-g0/Rgas/L[1])

        prss = (self.temperature/T[0])**(-g0/Rgas/L[0])*layer[0]
        prss = prss + P1*np.exp(-g0*(gpa-z[1])/Rgas/T[1])*layer[1]
        prss = prss + P2*np.exp(-g0*(gpa-z[2])/Rgas/T[2])*layer[2]

        return prss*101325
    
np.set_printoptions(precision=10, suppress=True)
hs = np.arange(0, 80000, 50)

satmos = stdday(hs)
catmos = coldday(hs)
hatmos = hotday(hs)
tatmos = tropday(hs)

plt.plot(satmos.pressure/101325, hs, '-o', label='S')
plt.plot(catmos.pressure/101325, hs, '-o', label='C')
plt.plot(hatmos.pressure/101325, hs, '-o', label='Hot')
plt.plot(tatmos.pressure/101325, hs, '-o', label='Trp')
plt.legend()
plt.show(block=False)

plt.figure()
plt.plot(satmos.temperature, hs, '-o', label='S')
plt.plot(catmos.temperature, hs, '-o', label='C')
plt.plot(hatmos.temperature, hs, '-o', label='Hot')
plt.plot(tatmos.temperature, hs, '-o', label='Trp')
plt.legend()
plt.show()
