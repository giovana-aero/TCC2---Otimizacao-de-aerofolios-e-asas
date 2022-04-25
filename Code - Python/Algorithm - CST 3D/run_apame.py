import numpy as np
import os
from math import isnan

def run_apame(pop,dat):
    
    # Esta função escreve o arquivo de entrada do APAME, faz a simulação da asa
    # e retorna os coeficientes aerodinâmicos relevantes
    
    # Dados da simulação
    cases = dat.cases
    v_ref = dat.v_ref
    rho = dat.rho
    p_atm = dat.p_atm
    aoa = dat.aoa
    mach = dat.mach
    # Dados da geometria
    b = pop.b
    mac = pop.mac
    NODE = pop.NODE
    PANEL = pop.PANEL
    sec_N = pop.sec_N
    S = pop.S
    
    
    aero = np.zeros((cases,4))
    for i in range(cases):
        
        # Imprimir o arquivo de entrada ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fid = open('apame_input.inp','w')
        fid.write('APAME input file\nVERSION 3.1\n\n\n')
        # Parâmetros do escoamento
        fid.write('# FLOW PARAMETERS\n\n')
        fid.write('# Airspeed [m/s]\n')
        fid.write('AIRSPEED ' + str(v_ref[i]) + '\n')
        fid.write('# Air density [kg/m^3\n')
        fid.write('DENSITY ' + str(rho[i]) + '\n')
        fid.write('# Atmospheric pressure [Pa]\n')
        fid.write('PRESSURE ' + str(p_atm[i]) + '\n')
        fid.write('# Prandtl-glauert corretion:\n#  0-no correction\n#  *-Mach number\n')
        fid.write('MACH ' + str(mach[i]) + '\n')
        fid.write('# Number of cases\n# Angle of attack [degrees]\n# sideslip angle [degrees]\n')
        fid.write('CASE_NUM 1\n')
        fid.write(str(aoa[i]) + '\n0\n\n\n')
        # Valores de referência
        fid.write('# REFERENCE VALUES\n\n')
        fid.write('# Wingspan [m]\n')
        fid.write('WINGSPAN ' + str(b) + '\n')
        fid.write('# Mean aerodynamic chord [m]\n')
        fid.write('MAC ' + str(mac) + '\n')
        fid.write('# Wing surface area [m^2]\n')
        fid.write('SURFACE ' + str(S) + '\n')
        fid.write('# Reference point [m]\n')
        fid.write('ORIGIN *\n0 0 0\n\n\n')
        # Parâmetros do solucionador
        fid.write('# SOLVER PARAMETERS\n\n')
        fid.write('# Singularity method:\n#  0-constant source/doublet\n#  1-constant doublet\n')
        fid.write('METHOD 0\n')
        fid.write('# Error\n')
        fid.write('ERROR 0.0000001\n')
        fid.write('# Collocation point depth\n')
        fid.write('COLLDIST 0.0000001\n')
        fid.write('# "Far field" coefficient\n')
        fid.write('FARFIELD 5\n')
        fid.write('# Collocation point calculation:\n#  0-approximative\n#  1-accurate\n')
        fid.write('COLLCALC 0\n')
        fid.write('# Interpolation method/order for velocity calculations:\n#  0-nodal\n#  1-first\n#  2-second\n')
        fid.write('VELORDER 1\n\n\n')
        # Resultados
        fid.write('# RESULT REQUESTS\n#  0-no\n# 1-yes\n')
        fid.write('RESULTS 1\n')
        fid.write('#  1 - coefiicients\n')
        fid.write('RES_COEF 1\n')
        fid.write('#  2 - forces\n')
        fid.write('RES_FORC 0\n')
        fid.write('#  3 - geometry\n')
        fid.write('RES_GEOM 0\n')
        fid.write('#  4 - velocity\n')
        fid.write('RES_VELO 0\n')
        fid.write('#  5 - pressure\n')
        fid.write('RES_PRES 0\n')
        fid.write('#  6 - center points\n')
        fid.write('RES_CENT 0\n')
        fid.write('#  7 - dipole values\n')
        fid.write('RES_DOUB 0\n')
        fid.write('#  8 - source values\n')
        fid.write('RES_SORC 0\n')
        fid.write('#  9 - velocity components\n')
        fid.write('RES_VELC 0\n')
        fid.write('# 10 - mesh characteristics\n')
        fid.write('RES_MESH 0\n')
        fid.write('# 11 - static pressure\n')
        fid.write('RES_STAT 0\n')
        fid.write('# 12 - dynamic pressure\n')
        fid.write('RES_DYNA 0\n')
        fid.write('# 13 - manometer pressure\n')
        fid.write('RES_MANO 0\n\n\n')
        # Geometria
        fid.write('# GEOMETRY\n\n')
        fid.write('NODES ' + str(NODE.shape[0]) + '\n')
        np.savetxt(fid, NODE, delimiter=' ', newline='\n', header='', footer='', comments='# ')
        fid.write('\nPANELS ' + str(PANEL.shape[0]) + '\n')
        np.savetxt(fid, PANEL[0:-sec_N+1,:], fmt='%d', delimiter=' ', newline='\n', header='', footer='', comments='# ')
        np.savetxt(fid, PANEL[-sec_N+1:,0:7], fmt='%d', delimiter=' ', newline='\n', header='', footer='', comments='# ')
        fid.close()
    
        # Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        # Nota: para fazer a simulação, é necessário abrir o executável e inserir o nome
        # do arquivo a ser analisado sem a extensão .inp. Devido à natureza desse processo,
        # é necessário o uso de um arquivo permanente contendo o nome do arquivo de input
        # do apame. A omissão da exetnsão do arquivo permanente foi realizada com o 
        # intuito de identificar melhor esse arquivo.
        os.system('apame_win64.exe < apame_cmd_input')
    
        # Apagar arquivos desnecessários
        os.remove('fort.2')
        # os.remove('apame_input.res')
    
        # Ler os resultados ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # fidpolar = open('apame_input.log')
        # dataBuffer1 = textscan(fidpolar,'#f #f #f #f #f #f','HeaderLines',51,... # CX, CY, CZ, CL, CM, CN (body reference frame)
        #     'CollectOutput',1,...
        #     'Delimiter',''); fclose(fidpolar);
        # fidpolar = fopen('apame_input.log');
        # dataBuffer2 = textscan(fidpolar,'#f #f #f','HeaderLines',55,... # CD, CK, CL (aerodynamic reference frame)
        #     'CollectOutput',1,...
        #     'Delimiter',''); fclose(fidpolar);
        
        dataBuffer1 = np.loadtxt('apame_input.log', skiprows=51, max_rows=1) 
        dataBuffer2 = np.loadtxt('apame_input.log', skiprows=55, max_rows=1) 
        
        CL = dataBuffer2[2]
        CD = dataBuffer2[0]
        CM = dataBuffer1[4]
            
        if isnan(CL) or CD <= 0:
            aero = 'n'; break
        else:
            aero[i,:] = [CL,CD,CL/CD,CM]
        
        
    
    
    # #delete('apame_input.log');
    
    return aero
    # return dataBuffer1,dataBuffer2