import os
import numpy as np

def run_xfoil_naca4_TCC2(naca_num,dat):
    
    reynolds = dat.reynolds
    aoa = dat.aoa
    iter_sim = dat.iter_sim
    cases = dat.cases

    # a_polar = 'polar.txt'

    # Apagar arquivos caso existam
    for i in range(cases):
        a_polar = 'polar' + str(i+1) + '.txt'
        if os.path.exists(a_polar):
            os.remove(a_polar)

    # Criar arquivo de input do XFOIL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    fid = open('xfoil_input.txt',"w") 

    # Mudar uma opção dos gráficos do XFOIL (desativar a aparição 
    # da janela com o desenho da simulação)
    fid.write('PLOP\nG\n\n')
    # fid.write('')

    # Nome do aerofólio
    fid.write('NACA ' + naca_num + '\n')

    # Simulações
    fid.write('OPER\n')
    
    if cases == 1: # Apenas uma condição de voo
        
        fid.write('VISC ' + str(reynolds[0]) + '\n') # Aplicar o modo viscoso e determinar número de Reynolds
        fid.write('PACC\n')              # Estabelecer arquivo de output
        fid.write('polar' + str(1) + '.txt\n\n')
        fid.write('ITER '+ str(iter_sim[0]) + '\n') # Mudar o número de iterações
        fid.write('ALFA ' + str(aoa[0]) + '\n')        # Estabelecer ângulo de ataque
    
    
    else: # Duas ou mais condições de voo
        
        for i in range(cases):
            
            if i == 0:
                fid.write('VISC ' + str(reynolds[i]) + '\n') # Aplicar o modo viscoso e determinar número de Reynolds
            elif i != 0 and reynolds[i] != reynolds[i-1]:
                fid.write('RE' + str(reynolds[i]) + '\n')    
            fid.write('PACC\n')              # Estabelecer arquivo de output
            fid.write('polar' + str(i+1) + '.txt\n\n')
            if i == 0 or i != 0 and iter_sim[i] != iter_sim[i-1]:
                fid.write('ITER' + str(iter_sim[i]) + '\n')
            fid.write('ALFA ' + str(aoa[i]) + '\n')        # Estabelecer ângulo de ataque
            fid.write('PACC\n') # fechar a polar para começar a próxima simulação

    # Fechar arquivo
    fid.write('\nQUIT\n')
    fid.close()

    # Executar XFOIL com o arquivo de entrada
    os.system("xfoil.exe < xfoil_input.txt") 

    # Ler o arquivo de saída ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aero = np.zeros((cases,4))
    for i in range(cases):
    
        # Contar o número de linhas
        nl = sum(1 for line in open('polar' + str(i+1) + '.txt'))
    
        if nl == 13:
            dataBufferPol = np.loadtxt('polar' + str(i+1) + '.txt', skiprows=12) 
    
            # Valores dos coeficientes
            CL = dataBufferPol[1]
            CD = dataBufferPol[2]
            CM = dataBufferPol[4]
    
            aero[i,:] = [CL,CD,CL/CD,CM]
    
        else:
            aero = 'n'
            break


    return aero