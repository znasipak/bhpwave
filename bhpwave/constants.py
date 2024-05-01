Modot_MKS = 1.98841e+30 # kg
GM_MKS = 1.32712440041279419e+20 # m^3/s^2
c_MKS = 299792458. # m/s
pc_MKS = 3.0856775814913674e+16 # m
yr_MKS = 31558149.763545603 # s (sidereal year)

Modot_GC1_to_S = GM_MKS/c_MKS**3
Modot_GC1_to_M = GM_MKS/c_MKS**2
Modot_GC1_to_PC = Modot_GC1_to_M/pc_MKS

A_MAX = 0.9999
OMEGA_MIN = 2.e-3