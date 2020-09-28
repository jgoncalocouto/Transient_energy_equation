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
geometry['H']=2*(2/3)
# Ambient Conditions
external={}
external['emissivity']=1.0
external['v_Wind']=15
T_amb=np.array([-40,-30,-25,-20,-15,-10,-5,0])

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
        external['T_film'] = (internal['T_in'] + external['T_amb']) * 0.5
        external['mu'] = thermal_properties('mu', 'air', external['T_film'])  # [m^2/(s)]
        external['rho'] = thermal_properties('rho', 'air', external['T_film'])  # [kg/(m^3)]
        external['k'] = thermal_properties('k', 'air', external['T_film']) # [W/m*°K]
        external['Pr'] = thermal_properties('Pr', 'air', external['T_film'])

            #Front Face
        external_facefront=external.copy()
        external_facefront['L_characteristic']=geometry['H']
        external_facefront['Re']=(external_facefront['rho']*external_facefront['v_Wind']*external_facefront['L_characteristic'])/external_facefront['mu']
        external_facefront=nusselt_external_frontplate(external_facefront)
        external_facefront['h'] = (external_facefront['Nu'] * external_facefront['k']) / (external_facefront['L_characteristic'])
        external_facefront['A_sup']=geometry['H']*geometry['L']

            #Laterals
        external_laterals = external.copy()
        external_laterals['L_characteristic'] = geometry['W']
        external_laterals['Re']=(external_laterals['rho']*external_laterals['v_Wind']*external_laterals['L_characteristic'])/external_laterals['mu']
        external_laterals=nusselt_external_flatplate(external_laterals)
        external_laterals['h'] = (external_laterals['Nu'] * external_laterals['k']) / (external_laterals['L_characteristic'])
        external_laterals['A_sup'] = 2*geometry['H'] * geometry['W']

            #Top Face
        external_top = external_laterals.copy()
        external_top['A_sup'] = geometry['L'] * geometry['W']

           #Back Face
        external_faceback=external.copy()
        external_faceback['L_characteristic']=geometry['H']
        external_faceback['Re']=(external_faceback['rho']*external_faceback['v_Wind']*external_faceback['L_characteristic'])/external_faceback['mu']
        external_faceback=nusselt_external_backplate(external_faceback)
        external_faceback['h'] = (external_faceback['Nu'] * external_faceback['k']) / (external_faceback['L_characteristic'])
        external_faceback['A_sup'] = geometry['H'] * geometry['L']

        external['A_sup']=(external_facefront['A_sup']+external_laterals['A_sup']+external_top['A_sup']+external_faceback['A_sup'])
        external['h_ext']=((external_facefront['h']*external_facefront['A_sup'])+(external_laterals['h']*external_laterals['A_sup'])+(external_top['h']*external_top['A_sup'])+(external_faceback['h']*external_faceback['A_sup']))/(external['A_sup'])

        # Calculate ventilation needs:
        ventilation = {}
        external['Q_rademit'] = sigma * external['A_sup'] * external['emissivity'] * (
                    ((internal['T_in'] + 273.15) ** 4) - ((external['T_amb'] + 273.15) ** 4))
        external['Q_conv'] = external['h_ext'] * (external['A_sup']) * (internal['T_in'] - external['T_amb'])

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

