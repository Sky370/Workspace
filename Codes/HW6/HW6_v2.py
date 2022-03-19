import numpy as np

# data
alp = 35
az1 = 310
az2 = 302.5
k = 5

# solution
a = np.sin(np.deg2rad(alp))*np.cos(np.deg2rad(az2-az1))
b = np.cos(np.deg2rad(alp))
phi = np.arctan(a/b)
phi_n = phi + np.arccos(np.cos(np.deg2rad(k))/np.sqrt(a**2+b**2))
phi_n = np.rad2deg(phi_n)

gamma = np.arccos((np.cos(np.deg2rad(alp))*np.cos(np.deg2rad(k))-np.cos(np.deg2rad(phi_n)))/(np.sin(np.deg2rad(k))*np.sin(np.deg2rad(alp))))
gamma = np.rad2deg(gamma)
