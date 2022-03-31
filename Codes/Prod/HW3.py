
import os
import pandas as pd
import numpy as np
from mpmath import *
import matplotlib.pyplot as plt


THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file1 = os.path.join(THIS_FOLDER, 'pressure.csv')
my_file2 = os.path.join(THIS_FOLDER, 'flowrate.csv')

# %%%%%%% Well testing HW5 %%%%%%%%%%%
k=50 #md
phi=0.25#
ct=1.5e-5 #psi-1
mu=1.2 #cp
B=1.25 #rb/stb
L=2000 #ft
A=2e4 #ft^2
Pwf=500  #psi
P0=3500 #psi
 
alpha=5.615
beta= float(1.127e-3)
yita=k/(phi*ct*mu)
 
time_rd = [0.0001, 0.01, 1, 100, 1000, 10000] #days need to trans to hours
x_rd = [0, 5, 10, 25, 50, 75, 100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000] #ft
 
# temp=sqrt(s/(alpha*beta*yita))
PressureAlongX = np.zeros((len(time_rd), len(x_rd)),dtype='float64')
FlowrateAlongX = np.zeros((len(time_rd), len(x_rd)),dtype='float64')
temp=beta*A*k/(B*mu)
for i in range(len(time_rd)):
    for j in range(len(x_rd)):
        fp = lambda s: P0/s-(P0-Pwf)/s*(cosh(sqrt(s/(alpha*beta*yita))*(x_rd[j]-L))/cosh(sqrt(s/(alpha*beta*yita))*L))
        PressureAlongX[i][j] = invertlaplace(fp, time_rd[i]/24, method='stehfest', degree = 12)
        fp = lambda s: -temp*(P0-Pwf)/s*sqrt(s/(alpha*beta*yita)*sinh(sqrt(s/(alpha*beta*yita))*(L-x_rd[j])))/cosh(sqrt(s/(alpha*beta*yita))*L)
        FlowrateAlongX[i][j] = invertlaplace(fp, time_rd[i]/24, method='talbot')
QFile =os.path.join(THIS_FOLDER, 'test.txt')
file = open(QFile, "w+")
pd.DataFrame(data=PressureAlongX).to_csv(my_file1)
pd.DataFrame(data=FlowrateAlongX).to_csv(my_file2)
content = str(PressureAlongX)
file.write(content)

a_file = open("test.txt", "w")
for row in PressureAlongX:
    np.savetxt(file, row)
print(PressureAlongX)
for i in range(len(time_rd)): 
    plt.plot(PressureAlongX[i][:])

plt.plot(x_rd, PressureAlongX[1][:])
plt.plot(x_rd, FlowrateAlongX[1][:])
