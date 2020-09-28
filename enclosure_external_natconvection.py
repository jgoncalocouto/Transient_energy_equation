#region Import Modules
from fluid_properties import *

import numpy as np
import pandas as pd
import math
from pyXSteam.XSteam import XSteam
import matplotlib.pyplot as plt
import pprint

#endregion

#region Inputs:

# Geometry
geometry={}
geometry['L']=0.69
geometry['W']=0.4
geometry['H']=2
# Ambient Conditions
external={}
external['emissivity']=0.94
T_amb=np.array([-50,-40,-30,-20,-10,0])
T_amb=np.array([-30])

# Internal Volume
internal={}
internal['T_in']=1

#endregion

#region Calculations

fig = plt.figure()
ax1 = fig.add_subplot(111)
Q_conv=np.zeros(len(T_amb))
Q_rademit=np.zeros(len(T_amb))
Q_balance=np.zeros(len(T_amb))

i=0
j=0
for T_ambi in T_amb:

        # Ambient Conditions
        external['A_sup'] = 2 * ( geometry['L'] * geometry['W'] + geometry['L'] * geometry['H'] + geometry['W'] * geometry['H'])  # [m^2]
        external['T_amb'] = T_ambi
        sigma = 5.67 * 10 ** (-8)

        # Internal Volume
        internal['V_total'] = geometry['H'] * geometry['L'] * geometry['W']

        # Conduction

        # Calculate convection coefficient
        g = 9.80665  # [m/s^2]
        external['T_film'] = (internal['T_in'] + external['T_amb']) * 0.5
        external['beta'] = 1 / (273.15 + external['T_film'])
        external['vu'] = thermal_properties('vu', 'air', external['T_film'])  # [m^2/(s)]
        external['k'] = thermal_properties('k', 'air', external['T_film'])  # [W/(m*ºK)]
        external['Pr'] = thermal_properties('Pr', 'air', external['T_film'])
        external['Gr'] = (g * external['beta'] * (geometry['H'] ** 3) * (internal['T_in'] - external['T_amb'])) / (
                    external['vu'] ** 2)
        external = nusselt_external_free(external)
        external['h_ext'] = (external['Nu'] * external['k']) / (geometry['H'])

        # Calculate ventilation needs:
        ventilation = {}
        external['Q_rademit'] = sigma * external['A_sup'] * external['emissivity'] * (
                    ((internal['T_in'] + 273.15) ** 4) - ((external['T_amb'] + 273.15) ** 4))
        external['Q_conv'] = external['h_ext'] * (external['A_sup']-2*geometry['W']*geometry['L']) * (internal['T_in'] - external['T_amb'])

        Q_balance[i]=(external['Q_conv']+external['Q_rademit'])
        Q_conv[i] = external['Q_conv']
        Q_rademit[i] = external['Q_rademit']
        i+=1
ax1.scatter(T_amb, Q_balance, label=str('Internal Temperature = '+ str(internal['T_in']) + ' [℃]'))
i=0

#endregion

#region Plots

# Plot1
ax1.set_xlabel('Ambient Temperature - T_{amb} - [ºC]')
ax1.set_xticks(T_amb)
ax1.set_ylabel('Thermal Power transferred to the exterior - [kW]')
ax1.set_title('Thermal Power Required to maintain internal temperature at T= '+str(internal['T_in'])+' [ºC]')
ax1.legend()
fig.show()

#endregion

