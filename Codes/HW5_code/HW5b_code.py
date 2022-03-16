from cmath import pi
from numpy import sqrt

# Given data
ROP = 50
N = 50
Angle = 65
D_o = 6.5
D_in = 4.5 
rho_f = 10
mu = 20
d_c = 0.1

# Calculation
v_c = 1/((1-(D_in/D_o)**2)*(0.64+18.16/ROP))
v_cs = 0.00516*mu+3.006

# C calculation
C_inc = 0.0342*Angle-0.000233*Angle**2-0.213
C_size = -1.04*d_c+1.286
C_mw = 1 - 0.0333*(rho_f-8.7)
C_geo_d = 0.277*(D_o-D_in)+0.2696
C_geo_mu = -0.00205841*mu+1.01493
C_geo_angle = 1.15667*10**-5*Angle**3-0.0026645*Angle**2+0.200266*Angle-3.9109
IMP = 0.555189 + 1.662083*N - 0.0224712*N*N + N*N*N*1.44*10**(-4)
C_r = (100-IMP)/100

v_s = v_cs*C_inc*C_size*C_mw*C_geo_d*C_geo_mu*C_geo_angle*C_r
v_f = v_s + v_c

# With rotation
Q = v_f*pi*((D_o/12)**2-(D_in/12)**2)/4*7.48*60

# With no rotation
C_r = 1
v_s = v_cs*C_inc*C_size*C_mw*C_geo_d*C_geo_mu*C_geo_angle*C_r
v_f = v_s + v_c
Q1 = v_f*pi*((D_o/12)**2-(D_in/12)**2)/4*7.48*60

print("The minimum flow rate that required to prevent bed development is: \nWith rotation - {q} [gpm] \nWith no rotation - {r} [gpm]".format(q = round(Q, 5), r = round(Q1, 5)))
