import numpy as np

# All are in SI unit
rho_l = 1000.5
rho_g = 1.225 
rho1 = 17.72
rho2 = 20.47
q_l = 0.00221 
q_g = 0.8495 
mw = 0.029
D_o = 0.2032
D_in = 0.127
D_h = D_o-D_in
H = 3048
dT_dL = 0.0273
T_sur = 294.26
P1 = 482633
k_ve = 0 #???????????
n_ve = 0.351
ty_ve = 0
R = 8.314

# Calculation
P_del = 68947.6
m_l = rho_l*q_l
m_g = rho_g*q_g
m_t = m_g+m_l
w_g = m_g/m_t
P2 = P_del+P1
T_diff = 10
T_int = T_sur
counter = 0

while T_diff > 0:
    a = w_g*R*T_int/mw
    b = (1-w_g)/rho_l
    rho = (b*(P_del)+a*np.log((a+b*P1)/(a+b*P2)))/b**2/(P_del)
    A = np.pi*(D_o**2 - D_in**2)/4
    v1 = m_t/(A*rho1)
    v2 = m_t/(A*rho2)
    v_f = m_t/(A*rho)
    eps = rho_l/rho
    k_ve = 0.2*(rho_l/rho)**(1-n_ve)
    D_e = D_h*3*n_ve/(2*n_ve+1)
    g_e = 12*v_f/D_e
    t_w = k_ve*(g_e)**n_ve
    L_del = (P_del)*D_h/(4*t_w+rho*9.81*D_h)
    T_new = T_sur+L_del*dT_dL
    T_diff = T_new - T_int
    T_int = T_new
    counter += 1

    print("Iteration {}, The section length L is: {}, ".format(counter, L_del), "Convergence is: {}".format(T_diff))
