'''
Created on Thu Feb 27 18:11:17 2022
Author: elh1873
'''

import pandas as pd
import numpy as np


# ----- >>>>>>
# Problem 2 (a)
# We will try to use pandas in order to use the .xlsx file for our project

df = pd.read_excel('Data_HW2.xlsx')
# print (df)

mol = df["Moles, (mol)"].tolist()
bp = df["B.P.\n(˚F)"].tolist()
mn = df["M\n(lbm/lbmole)"].tolist()
pc= df["Pc (psia)"].tolist()
tc = df["Tc (˚F)"].tolist()
V_c_raw = df['V_c_raw'].tolist()
P_ch = df['P_ch'].tolist()
w = df["ω"].tolist()
tf = df["Tr"].tolist()
n = len(mol)                                 # Number of components
P = 1000                                     # psia
R = 10.73157                                 # ˚R
T_f = 90                                     # ˚F
T_r = T_f + 459.67                           # ˚R

# Wilson's correlation
K = [(pc[i]/P)*np.exp(5.37*(1+w[i])*(1-tf[i]/T_r)) for i in range(n)]

# Peng Robinson correlation
Ca = 0.45724
Cb = 0.07780
v_f = float(input("Please enter your v/f value: "))

Nfg = 1
Cvg = 1
counter = 1

while Nfg > 10**-6 or Cvg > 10**-4:
    Z = [mol[i]/100 for i in range(n)]
    rob_t = [T_r/tf[i] for i in range(n)]
    rob_s = [0.37464+1.5422*(w[i])-0.26992*(w[i]**2) for i in range(n)]
    alpha = [(1+rob_s[i]*(1-np.sqrt(rob_t[i])))**2 for i in range(n)]
    x_i = [Z[i]/(1+v_f*(K[i]-1)) for i in range(n)]
    y_i = [Z[i]*K[i]/(1+v_f*(K[i]-1)) for i in range(n)]

    # Denominators
    dnomin1 = [x_i[i]*np.sqrt(alpha[i])*tf[i]/np.sqrt(pc[i])for i in range(n)]
    dnomin2 = [y_i[i]*np.sqrt(alpha[i])*tf[i]/np.sqrt(pc[i])for i in range(n)]
    dnomin3 = [x_i[i]*tf[i]/pc[i] for i in range(n)]
    dnomin4 = [y_i[i]*tf[i]/pc[i] for i in range(n)]

    # Coefficients
    ai_a_L = [np.sqrt(alpha[i])*tf[i]/np.sqrt(pc[i])/sum(dnomin1) for i in range(n)]    # ai/aL
    ai_a_V = [np.sqrt(alpha[i])*tf[i]/np.sqrt(pc[i])/sum(dnomin2) for i in range(n)]    # ai/aV
    bi_b_L = [tf[i]/pc[i]/sum(dnomin3) for i in range(n)]                               # bi/bL
    bi_b_V = [tf[i]/pc[i]/sum(dnomin4) for i in range(n)]                               # bi/bV
    A_L = Ca*P*(sum(dnomin1)/T_r)**2
    A_V = Ca*P*(sum(dnomin2)/T_r)**2
    B_L = Cb*P*sum(dnomin3)/T_r
    B_V = Cb*P*sum(dnomin4)/T_r

    mol_frac = [Z[i]*(K[i]-1)/(v_f*(K[i]-1)+1) for i in range(n)]
    mol_frac2 = [Z[i]*((K[i]-1)**2)/((v_f*(K[i]-1)+1)**2) for i in range(n)]
    v_f_new = v_f-(sum(mol_frac)/-sum(mol_frac2))

    # Determining Z_L and Z_V
    # Robinson Cubic Equation
    eq_ZL = [1, -(1-B_L), A_L-2*B_L-3*B_L**2, -(A_L*B_L-B_L**2-B_L**3)]
    eq_ZV = [1, -(1-B_V), A_V-2*B_V-3*B_V**2, -(A_V*B_V-B_V**2-B_V**3)]

    # Robinson roots
    if len(np.roots(eq_ZL)) < 2:                            # Single Phase
        root_s = np.roots(eq_ZL).real
    else:
        for i in range(len(np.roots(eq_ZL))):               # Multi phase
            if abs(np.roots(eq_ZL)[i].imag) < 10**-5:
                Z_L = np.roots(eq_ZL)[i].real
        for j in range(len(np.roots(eq_ZV))):
            if abs(np.roots(eq_ZV)[j].imag) < 10**-5:
                Z_V = np.roots(eq_ZV)[j].real

    # Fugacity
    fug_L = [np.exp(bi_b_L[i]*(Z_L-1)-np.log(Z_L-B_L)-A_L/B_L*np.log(1+B_L/Z_L)*(2*ai_a_L[i]-bi_b_L[i])) for i in range(n)]
    fug_V = [np.exp(bi_b_V[i]*(Z_V-1)-np.log(Z_V-B_V)-A_V/B_V*np.log(1+B_V/Z_V)*(2*ai_a_V[i]-bi_b_V[i])) for i in range(n)]

    # K from iterations
    K_new = [fug_L[i]/fug_V[i] for i in range(n)]

    # Convergence
    Cvg = sum([(K[i]/K_new[i]-1)**2 for i in range(n)])     # Convergence for K
    Nfg = v_f - v_f_new                                     # Convergence for v_f
    v_f = v_f_new
    K = K_new
    
    
    print('Iteration {},  Value of v/f : {},  Value of Convergence K: {}'.format(counter, str(v_f), str(Cvg)))
    counter += 1
    
    # End of Loop

print('\nThe final value of v/f is: {}\n'.format(str(v_f)))

# Phase envelop
Temps = list(range(200, 900, 50))
Denomin = [Z[i]*(1-tf[i]/Temps[i])/(pc[i]*np.exp(1+w[i])) for i in range(n)]

# ---------    Stage 2    -------------

# Problem 2 (b)

# Density calculation
Rho_L = sum([P*x_i[i]*mn[i]/(Z_L*R*T_r) for i in range(n)])
Rho_V = sum([P*y_i[i]*mn[i]/(Z_V*R*T_r) for i in range(n)])

# Liq and Gas viscosity
Tpc = [tf[i]*Z[i] for i in range(n)]    # Pseudo Critical Temperature
Ppc = [pc[i]*Z[i] for i in range(n)]    # Pseudo Critical Pressure
P_b = [Z[i]*pc[i]*np.exp(5.37*(1+w[i])*(1-tf[i]/T_r)) for i in range(n)]        # Bubble Pressure
eps = [5.4402*tf[i]**(1/6)/(np.sqrt(mn[i])*pc[i]**(2/3)) for i in range(n)]     # epsilon
myu = [34*10e-6*rob_t[i]**0.94/eps[i] if rob_t[i] <= 1.5 else 17.78*10e-6*(((4.58*rob_t[i])-1.67)**0.625)/eps[i] for i in range(n)]


Vc_7plus = [0.463e-3, 0.667e-3, 0.92e-3, 1.67e-3 ]  # For components after C7
M7_plus =  [143e-6, 166e-6, 230e-6, 409e-6]
x_i_c7= x_i[10:14]
Sg_C7plus = [M7_plus[i] / Vc_7plus[i] for i in range(len(Vc_7plus))]
Sg_C7 = sum([x_i_c7[i]*Sg_C7plus[i] for i in range(len(Vc_7plus))])
VC_7_plus = 21.573 + 0.015122*(sum(M7_plus))-27.656*Sg_C7+0.070615*(sum(M7_plus))*Sg_C7

Y_ol = sum([x_i[i]*np.sqrt(mn[i])*myu[i] for i in range(n)])/sum([x_i[i]*np.sqrt(mn[i]) for i in range(n)])
V_c = [V_c_raw[i]/mn[i] for i in range(n)]
ximiVci = [x_i[i]*mol[i]*V_c[i] for i in range(10)]

# Calculating viscosity
rho_reduced = Rho_L/sum(mol)*sum(ximiVci)+sum(x_i_c7)*sum(Vc_7plus)
esp_m = 5.4402*((sum(Tpc)**(1/6)))/(((sum(mol))**0.5)*((sum(Ppc))**(2/3)))
Visc = Y_ol + (esp_m**(-1))*((0.1023+0.023364*rho_reduced+0.058533*rho_reduced**2-0.040758*rho_reduced**3+0.0093724*rho_reduced**4)**4-0.0001)

# Surface tension
Coef1 = Rho_L/(62.4*sum(mol[10:]))
Coef2 = Rho_V/(62.4*sum(mol[:9]))
sig25 = [P_ch[i]*(Coef1*x_i[i]-Coef2*y_i[i]) for i in range(n)]

# Liveoil
lvoil = [mol[i]/V_c[i] for i in range(n)]

print('Gas density is: {} [lb/ft3]'.format(Rho_V))
print('Liquid density is: {} [lb/ft3]'.format(Rho_L))
print('Bubble pressure is: {} [psia]'.format(sum(P_b)))
print('Pseudo critical temperature is: {} [R]'.format(sum(Tpc)))
print('Pseudo critical pressure is: {} [Psia]'.format(sum(Ppc)))
print('Viscosity is: {} [cp]'.format(Visc))
print('Surface tension is: {} [dynes/cm]'.format(sum(sig25)))
print('Liveoil density is: {} [lb/ft3]'.format(sum(lvoil)/62.4))
