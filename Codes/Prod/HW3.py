'''
Created on Thu Mar 19 14:00:00 2022
Author: elh1873
'''

import pandas as pd
import numpy as np

df = pd.read_excel('Data_HW2.xlsx')

# Data from Excel
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

