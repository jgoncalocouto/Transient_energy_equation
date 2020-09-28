from fluid_properties import *

def external_natural_convection(T_in,T_amb,material,L_characteristic,A_ht):
    Q_conv=0
    material=material.lower()
    g = 9.80665  # [m/s^2]
    T_film = (T_in+T_amb) * 0.5
    beta = 1 / (273.15 + T_film)
    vu = thermal_properties('vu', material, T_film)  # [m^2/(s)]
    k = thermal_properties('k', material, T_film)  # [W/(m*ºK)]
    Pr = thermal_properties('Pr', 'air', T_film)
    Gr = (g * beta * (L_characteristic ** 3) * (T_in - T_amb)) / (vu ** 2)
    ht={}
    ht['Gr'] = Gr
    ht['Pr']=Pr
    ht = nusselt_external_free(ht)
    Nu= ht['Nu']
    h_ext = (Nu * k) / (L_characteristic)
    Q_conv = h_ext * (A_ht) * (T_amb-T_in)
    return Q_conv

def external_forced_convection(T_in,T_amb,material,V_wind=0.5,L=0.4,W=0.7,H=2):
    Q_conv=0
    external={}
    geometry={}
    geometry['H']=H; geometry['W']=W; geometry['L']=L;
    external['v_Wind'] = V_wind
    external['T_film']=(T_in + T_amb) * 0.5
    external['mu'] = thermal_properties('mu', material, external['T_film'])  # [m^2/(s)]
    external['rho'] = thermal_properties('rho', material, external['T_film'])  # [kg/(m^3)]
    external['k'] = thermal_properties('k', material, external['T_film'])  # [W/m*°K]
    external['Pr'] = thermal_properties('Pr', material, external['T_film'])

    # Front Face
    external_facefront = external.copy()
    external_facefront['L_characteristic'] = geometry['H']
    external_facefront['Re'] = (external_facefront['rho'] * external_facefront['v_Wind'] * external_facefront[
        'L_characteristic']) / external_facefront['mu']
    external_facefront = nusselt_external_frontplate(external_facefront)
    external_facefront['h'] = (external_facefront['Nu'] * external_facefront['k']) / (
        external_facefront['L_characteristic'])
    external_facefront['A_sup'] = geometry['H'] * geometry['L']

    # Laterals
    external_laterals = external.copy()
    external_laterals['L_characteristic'] = geometry['W']
    external_laterals['Re'] = (external_laterals['rho'] * external_laterals['v_Wind'] * external_laterals[
        'L_characteristic']) / external_laterals['mu']
    external_laterals = nusselt_external_flatplate(external_laterals)
    external_laterals['h'] = (external_laterals['Nu'] * external_laterals['k']) / (
        external_laterals['L_characteristic'])
    external_laterals['A_sup'] = 2 * geometry['H'] * geometry['W']

    # Top Face
    external_top = external_laterals.copy()
    external_top['A_sup'] = geometry['L'] * geometry['W']

    # Back Face
    external_faceback = external.copy()
    external_faceback['L_characteristic'] = geometry['H']
    external_faceback['Re'] = (external_faceback['rho'] * external_faceback['v_Wind'] * external_faceback[
        'L_characteristic']) / external_faceback['mu']
    external_faceback = nusselt_external_backplate(external_faceback)
    external_faceback['h'] = (external_faceback['Nu'] * external_faceback['k']) / (
    external_faceback['L_characteristic'])
    external_faceback['A_sup'] = geometry['H'] * geometry['L']

    external['A_sup'] = (
                external_facefront['A_sup'] + external_laterals['A_sup'] + external_top['A_sup'] + external_faceback[
            'A_sup'])
    external['h_ext'] = ((external_facefront['h'] * external_facefront['A_sup']) + (
                external_laterals['h'] * external_laterals['A_sup']) + (external_top['h'] * external_top['A_sup']) + (
                                     external_faceback['h'] * external_faceback['A_sup'])) / (external['A_sup'])

    external['Q_conv'] = external['h_ext'] * (external['A_sup']) * (T_amb-T_in)
    Q_conv=external['Q_conv']
    return Q_conv

def heating_resistance(T,Q,T_on=0,T_off=10):

    if T<=T_on:
        Q_heat=Q
    elif T>=T_off:
        Q_heat=0
    else:
        Q_heat=0


    return Q_heat