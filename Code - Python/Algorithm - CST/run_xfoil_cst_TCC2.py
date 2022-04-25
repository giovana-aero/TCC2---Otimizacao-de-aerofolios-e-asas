import os
import numpy as np

def run_xfoil_cst_TCC2(dat,op):
    
    if op == 1: # Montar o arquivo de input
    
        # Pegar as informações do código base
        reynolds = dat.reynolds
        aoa = dat.aoa
        iter_sim = dat.iter_sim
        numNodes = dat.numNodes
        cases = dat.cases

    # a_polar = 'polar.txt'

    # if os.path.exists('xfoil_input.txt'):
    #     os.remove('xfoil_input.txt')

        # Criar arquivo de input do XFOIL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        fid = open('xfoil_input.txt',"w") 

        # Mudar uma opção dos gráficos do XFOIL (desativar a aparição 
        # da janela com o desenho da simulação)
        fid.write('PLOP\nG\n\n')

        # Ler coordenadas
        fid.write('LOAD coordenadas.dat\n\n') # Nota: retirar a quebra de linha extra
                                              # caso o perfil tenha um nome

        # Alterar número de nós
        if numNodes != 0:
            fid.write('PPAR\n')
            fid.write('N ' + str(numNodes) + '\n')
            fid.write('\n\n')
    
        # Simulação
        fid.write('OPER\n')
        if cases == 1: # Apenas uma condição de voo
            fid.write('VISC ' + str(reynolds[0]) + '\n') # Aplicar o modo viscoso
            fid.write('PACC\n')              # Estabelecer arquivo de output
            fid.write('polar' + str(1) + '.txt\n\n')
            fid.write('ITER ' + str(iter_sim[0]) + '\n') # Mudar o número de iterações
            fid.write('ALFA ' + str(aoa[0]) + '\n')        # Estabelecer ângulo de ataque
        
        else: # Duas ou mais condições de voo
            for i in range(cases):
                
                if i == 0:
                    fid.write('VISC ' + str(reynolds[i]) + '\n') # Aplicar o modo viscoso
                elif i != 0 and reynolds[i] != reynolds[i-1]:
                    fid.write('RE ' + reynolds[i] + '\n')
                fid.write('PACC\n')              # Estabelecer arquivo de output
                fid.write('polar' + str(i+1) + '.txt\n\n')
                if i == 0 or i != 0 and iter_sim[i] != iter_sim[i-1]:
                    fid.write('ITER ' + str(iter_sim[i]) + '\n') # Mudar o número de iterações
                fid.write('ALFA ' + str(aoa[i]) + '\n')        # Estabelecer ângulo de ataque
                fid.write('PACC\n')

    
        # Fechar arquivo
        fid.write('\nQUIT\n')
        fid.close()



    else: # Simular os aerofólios
        cases = dat.cases
        
        # Apagar arquivos caso existam
        for i in range(cases):
            a_polar = 'polar' + str(i+1) + '.txt'
            if os.path.exists(a_polar):
                os.remove(a_polar)
        
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