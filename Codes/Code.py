import math
from turtle import color
import matplotlib.pyplot as plt
import numpy as np


# Given Data
A = [0.5, 2, 0.7, 4, 0.4, 1.5]
B = [7e-8, 9e-8, 3e-8, 9e-7, 1.8e-8, 8e-8]
C = [1.5, 1.2, 3.0, 1.5, 1.0, 1.5]
Pwf_IPR = []
Pwf_OPR = []
q_rateM = list(range(0, 30000000, 1000))

# Adjustments
P_inj = 3800**2  # Injection pressure, psia
P_res = 235**2  # Reservoir pressure, psia
injection_rate_step = 1  # Injection rate increase by 1 x MMscf/d at each step

# Iterating through every element in A
for i in range(len(A)):
    q_rate = q_rateM
    int_point = 0
    q0 = 37000000
    p = 0
    for value in q_rate:
        p = p + 1 
        if A[i] * P_inj < B[i]*value*value:
            q0 = value
            break
        y = math.sqrt((value/C[i])+P_res)
        Pwf_OPR.append(y)
        x = math.sqrt(A[i]*P_inj-B[i]*value*value)
        Pwf_IPR.append(x)
        
        if (abs(x-y) < 1):
            int_point = p
    # Plotting the data
    q_rate = q_rateM[0:len(Pwf_IPR)]
    plt.title("Gas Flow Rate vs P_wf ", fontdict=None, loc='center', pad=None)
    plt.xlabel("Gas Flow Rate [MMscf/d]")
    plt.ylabel("P_wf [psia]")
    plt.plot(q_rate, Pwf_IPR, '#20b2aa', label="IPR")
    plt.plot(q_rate, Pwf_OPR, '#2F5233', label="OPR")
    plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, loc="upper right")
    plt.xlim(0, q0)
    plt.ylim(0,)

    # Intersection points and slopes are being calculated
    arr_q_rate = np.array(q_rate)
    arr_Pwf_IPR = np.array(Pwf_IPR)
    arr_Pwf_OPR = np.array(Pwf_OPR)
    idx = np.argwhere(np.diff(np.sign(arr_Pwf_IPR-arr_Pwf_OPR))).flatten()
    plt.plot(arr_q_rate[idx], arr_Pwf_IPR[idx], 'ob')
    print("Flow Rate  at Operational point: , {}, [MMscf/d],".format(arr_q_rate[idx]/1000000)
          + "\nP_wr @IPR: , {}, [psia]".format(arr_Pwf_IPR[idx]), "\nP_wr @OPR: , {}, [psia]".format(arr_Pwf_OPR[idx]))

    IPR_slope = abs((Pwf_IPR[int_point+10]-Pwf_IPR[int_point-10]) /
                    (q_rate[int_point+10]-q_rate[int_point-10]))
    OPR_slope = abs((Pwf_OPR[int_point+10]-Pwf_OPR[int_point-10]) /
                    (q_rate[int_point+10]-q_rate[int_point-10]))
    if IPR_slope > OPR_slope:
        print("IPR Slope > OPR Slope. Nodal analysis is stable!" +
              "\n\t |m| > |M|\n\t {} > {} ".format(str(round(IPR_slope, 5)), str(round(OPR_slope, 5))))
    else:
        print("IPR Slope < OPR Slope. Nodal analysis is unstable!" +
              "\n\t |m| < |M|\n\t {} < {} ".format(str(round(IPR_slope, 5)), str(round(OPR_slope, 5))))
    
    plt.show()
    Pwf_IPR.clear()
    Pwf_OPR.clear()