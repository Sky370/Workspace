import numpy as np
import matplotlib.pyplot as plt

# Given data
MD = [0, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]
inc = [0, 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 30.2, 30.4, 30.3, 30.6, 31, 31.2, 30.7, 31.4, 30.6, 30.5]
az = [0, 0, 21.7, 26.5, 23.3, 20.3, 23.3, 23.9, 24.4, 23.4, 23.7, 23.3, 22.8, 35.8, 47.6, 59.3, 65.2, 78.4, 91.2, 108.7, 110.2, 110.0] 

inc_rad = [np.deg2rad(inc[i]) for i in range(len(inc))]
az_rad = [np.deg2rad(az[i]) for i in range(len(az))]
T = [(az_rad[i+1]-az_rad[i])/(MD[i+1]-MD[i]) for i in range(len(MD)-1)]        # Turn Rate
B = [(inc_rad[i+1]-inc_rad[i])/(MD[i+1]-MD[i]) for i in range(len(MD)-1)]      # Build Rate
H = [(az[i+1]-az[i])/((MD[i+1]-MD[i])*np.sin((inc_rad[i]+inc_rad[i+1])/2)) for i in range(len(MD)-1)]
H = [np.deg2rad(H[i]) for i in range(len(H))]                                  # Horizontal Turn Rate

# Average Angle Method
del_x = [np.sin((inc_rad[i]+inc_rad[i+1])/2)*np.cos((az_rad[i]+az_rad[i+1])/2)*(MD[i+1]-MD[i]) for i in range(len(MD)-1)]
del_y = [np.sin((inc_rad[i]+inc_rad[i+1])/2)*np.sin((az_rad[i]+az_rad[i+1])/2)*(MD[i+1]-MD[i]) for i in range(len(MD)-1)]
del_z = [np.cos((inc_rad[i]+inc_rad[i+1])/2)*(MD[i+1]-MD[i]) for i in range(len(MD)-1)]

# Constant Build and Horizontal Turn Rate Method
del_x1 = [(np.sin(az_rad[i+1])-np.sin(az_rad[i]))/H[i] for i in range(len(MD)-1)]
del_y1 = [(np.cos(az_rad[i])-np.cos(az_rad[i+1]))/H[i] for i in range(len(MD)-1)]
del_z1 = [(np.sin(inc_rad[i+1])-np.sin(inc_rad[i]))/B[i] for i in range(len(MD)-1)]

# Constant Build and Turn Rate Method
del_x2 = [(T[i]*(np.sin(inc_rad[i+1])*np.sin(az_rad[i+1])-np.sin(inc_rad[i])*np.sin(az_rad[i]))+B[i]*(np.cos(inc_rad[i+1])*np.cos(az_rad[i+1])-np.cos(inc_rad[i])*np.cos(az_rad[i])))/(T[i]**2-B[i]**2) for i in range(len(MD)-1)]
del_y2 = [(B[i]*(np.sin(az_rad[i+1])*np.cos(inc_rad[i+1])-np.sin(az_rad[i])*np.cos(inc_rad[i]))-T[i]*(np.cos(az_rad[i+1])*np.sin(inc_rad[i+1])-np.cos(az_rad[i])*np.sin(inc_rad[i])))/(T[i]**2-B[i]**2) for i in range(len(MD)-1)]
del_z2 = [(np.sin(inc_rad[i+1])-np.sin(inc_rad[i]))/B[i] for i in range(len(MD)-1)]

# Constant Dogleg Method
betta = [2*np.arcsin(np.sqrt(np.sin((inc_rad[i+1]-inc_rad[i])/2)**2+np.sin(inc_rad[i])*np.sin(inc_rad[i+1])*np.sin((az_rad[i+1]-az_rad[i])/2)**2)) for i in range(len(MD)-1)]
RF = [(MD[i+1]-MD[i])/betta[i]*np.tan(betta[i]/2) for i in range(len(MD)-1)]
del_x3 = [RF[i]*(np.sin(inc_rad[i])*np.cos(az_rad[i])+np.sin(inc_rad[i+1])*np.cos(az_rad[i+1])) for i in range(len(MD)-1)]
del_y3 = [RF[i]*(np.sin(inc_rad[i])*np.sin(az_rad[i])+np.sin(inc_rad[i+1])*np.sin(az_rad[i+1])) for i in range(len(MD)-1)]
del_z3 = [RF[i]*(np.cos(inc_rad[i])+np.cos(inc_rad[i+1])) for i in range(len(MD)-1)]

fig = plt.figure()
ax = plt.axes(projection='3d')
# Data for a three-dimensional line
zline = del_z
xline = del_x
yline = del_y
ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
zdata = del_z
xdata = del_x
ydata = del_y
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')

plt.show()
