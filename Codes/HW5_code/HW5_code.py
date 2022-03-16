from cmath import pi
from numpy import sqrt

# Given data
Cc = 0.02
ROP = 50
D_o = 8.5
D_in = 4.5
rho_f = 10
mu = 10
tao = 5
d_c = 0.1
rho_c = 2.65*8.33

# Calculation
v_s = float(input("Input your slip velocity: "))
v_diff = 1
counter = 0

while v_diff > 10**-4:
    v_f = v_s + 2.785*10**-4*D_o**2*ROP/(Cc*(D_o**2-D_in**2))
    mu_e = mu + 5*tao*d_c/v_f
    N_re = rho_f*v_s*d_c*928/mu_e
    if N_re <= 1:
        C_d = 40/N_re
    elif 1<N_re<1000:
        C_d = 22/sqrt(N_re)
    else:
        C_d = 1.5
    v_s_new = 1.89*sqrt(d_c/C_d*(rho_c-rho_f)/rho_f)
    v_diff = v_s - v_s_new
    v_s = v_s_new
    counter += 1
    print("Iteration {}, Slip velocity is: {}".format(counter, v_s))


Q = ROP*(D_o)**2/(1466.95*Cc*(1-v_s/v_f))

print("The required flow rate is: {} [gpm]".format(round(Q, 5)))