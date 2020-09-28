#region Import Modules

from fluid_properties import *
from auxiliary_functions import *
import numpy as np
import pandas as pd
import math
from scipy import interpolate

from pyXSteam.XSteam import XSteam
import matplotlib.pyplot as plt
import pprint

#endregion
case = {'material': ['air'], 'environment_conditions': ['case_1'], 'L': [0.4], 'W': [0.7], 'H': [2], 'epsilon': [0.94], 'H_rad': [0.6], 'L_rad':[0.3]}
df_case = pd.DataFrame(case)


case_nr=0
f_veq=0.7;
T_on=0
T_off=10
Q_resistance=1000;


environment_xls = pd.ExcelFile('environment_conditions.xlsx')
df_environment = environment_xls.parse(df_case['environment_conditions'].loc[case_nr])

material=df_case['material'].loc[case_nr]
H=df_case['H'].loc[case_nr]
L=df_case['L'].loc[case_nr]
W=df_case['W'].loc[case_nr]
L_rad=df_case['L_rad'].loc[case_nr]
H_rad=df_case['L_rad'].loc[case_nr]
epsilon=df_case['epsilon'].loc[case_nr]
V=L*W*H*f_veq

#row = [1.4, 1.01, 0.4,0.7,2]  #adding new row to dataframe
#df_case.loc[len(df_case)] = row

#inputs:
t_step=60*10; #[s]
t_total=24; #[h]

t=np.arange(0, t_total*3600,t_step)
T=np.zeros((math.ceil(t_total*3600/t_step),1))
rho=np.zeros((math.ceil(t_total*3600/t_step),1))
cp=np.zeros((math.ceil(t_total*3600/t_step),1))
Q_rad=np.zeros((math.ceil(t_total*3600/t_step),1))
Q_rademit=np.zeros((math.ceil(t_total*3600/t_step),1))
Q_loss=np.zeros((math.ceil(t_total*3600/t_step),1))
Q_diss=np.zeros((math.ceil(t_total*3600/t_step),1))
Q_heating=np.zeros((math.ceil(t_total*3600/t_step),1))
T_amb=np.zeros((math.ceil(t_total*3600/t_step),1))
I_rad=np.zeros((math.ceil(t_total*3600/t_step),1))
v_wind=np.zeros((math.ceil(t_total*3600/t_step),1))

fi_tamb=interpolate.interp1d(df_environment['t']*3600, df_environment['T_amb'])
T_amb[:,0]=fi_tamb(t)
fi_irad=interpolate.interp1d(df_environment['t']*3600, df_environment['I_rad'])
I_rad[:,0]=fi_irad(t)
fi_vwind=interpolate.interp1d(df_environment['t']*3600, df_environment['v_wind'])
v_wind[:,0]=fi_vwind(t)



T[0]=T_amb[0] #Cabinet is in thermal equilibrium with environment

A_ht = 2 * ( L * W + L * H + W * H)  #[m^2]
A_rad=H_rad*L_rad
L_characteristic=H

for i in range(1,len(t)):
    cp[i-1]=thermal_properties('cp', material, T[i-1])
    rho[i - 1] = thermal_properties('rho', material, T[i - 1])
    cp[i-1]=1005
    rho[i-1]=1.2
    Q_heating[i]=heating_resistance(T[i-1],Q_resistance,T_on,T_off)
    #Q_loss[i]=external_natural_convection(T[i-1], T_amb[i], material, L_characteristic, A_ht)
    Q_loss[i]=external_forced_convection(T[i-1],T_amb[i],material,v_wind[i],L,W,H)
    Q_rademit[i]=epsilon * (5.67 * 10 ** (-8)) * A_ht * ((T_amb[i]+273.15) ** 4 - (T[i - 1]+273.15) ** 4)
    #Q_rad[i] = I_rad[i] * A_rad

    T[i]=T[i-1]+((Q_rad[i]+Q_loss[i]+Q_rademit[i]+Q_diss[i]+Q_heating[i])*((t[i]-t[i-1])/(rho[i-1]*V*cp[i-1]*1000)))
    if i % 360 == 0:
        print('Iteration: ' + str(i) + '; T= ' + str(T[i]))

#region Plots

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
fig.suptitle('1-D Transient Model')
ax1.plot(t, T,'b-')
ax1.plot(t, T_amb,'k-')
ax2.plot(t, Q_heating,'r')


#endregion
