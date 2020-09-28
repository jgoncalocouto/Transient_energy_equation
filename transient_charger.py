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
external['emissivity']=0.8
T_amb=np.array([-30])

# Internal Volume
internal={}
internal['duct_percentage_volume']=0.2
internal['T_in']=50

#endregion

#region Calculations
fig = plt.figure()
ax1 = fig.add_subplot(111)

i=0
j=0
for T_ambi in T_amb:
    for ti in t:


        # Ambient Conditions
        external['A_sup'] = 2 * ( geometry['L'] * geometry['W'] + geometry['L'] * geometry['H'] + geometry['W'] * geometry['H'])  # [m^2]
        external['T_amb'] = T_ambi

        # Internal Volume
        internal['V_total'] = geometry['H'] * geometry['L'] * geometry['W']

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
        ventilation['cp'] = thermal_properties('cp', 'air', internal['T_in']) * 1000  # [J/(kg*ºC)]
        ventilation['rho'] = thermal_properties('rho', 'air', internal['T_in'])  # [kg/(m^3)]
        ventilation['m_dot'] = ((power_eletronics['P_total'] + external['Q_radinc']) - (
                    external['Q_rademit'] + external['Q_conv'])) / (
                                           ventilation['cp'] * (internal['T_in'] - external['T_amb']))
        ventilation['V_dot'] = (ventilation['m_dot'] / ventilation['rho']) * 60  # [m^3/min]
        ventilation['rho_std'] = thermal_properties('rho', 'air', 15)  # [kg/(m^3)]
        ventilation['V_dot_std'] = (ventilation['m_dot'] / ventilation['rho_std']) * 3600  # [m^3/h]

        Vdot_std[i,j]=ventilation['V_dot_std']
        Q_conv[i,j]=external['Q_conv']
        Q_rademit[i,j]=external['Q_rademit']
        Q_balance[i,j]=external['Q_radinc']-(external['Q_conv']+external['Q_rademit'])

        i=i+1
    ax1.scatter(T_in, Vdot_std[:, j], label=str('Ambient Temperature = '+ str(T_ambi) + ' [℃]'))
    i=0
    j=j+1
#endregion
