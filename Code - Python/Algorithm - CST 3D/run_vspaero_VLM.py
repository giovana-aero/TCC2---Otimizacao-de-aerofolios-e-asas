import os
import numpy as np
# import csv
# import matplotlib.pyplot as plt
import multiprocessing
from time import sleep

# Esta função monta a asa no formato do OpenVSP e faz a simulação no VSPAERO

def execute():
    os.system('vspaero -omp 4 ' + 'wing_geom_DegenGeom') 


def run_vspaero_VLM(pop,dat):
    
    wait = 6 # Tempo limite para as simulações funcionarem [s]
    tries = 3 # Número máximo de tentativas de um mesmo indivíduo
    
    
    # Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    aero = np.zeros((dat.cases,4))
    
    for P in range(dat.cases):

        # Dados da geometria da asa
        plan_type = pop.type
        b = pop.b
        ST = pop.S
        mac = pop.mac
        if plan_type == 1:
            b1 = pop.b1
            b2 = b - b1

        # Dados da simulação
        aoa = dat.aoa[P] # Ângulos de ataque
        Vinf = dat.v_ref[P] # Velocidade
        rho = dat.rho[P] # Densidade do ar
        Re_ref = dat.reynolds[P] # Número de Reynolds na corda de referência
        KTC = dat.karman_tsien[P] # Opção de correção de Karman_Tsien
        wake = dat.wake_iters[P] # Número de iteações na simulação
        
        # # Mais alguns dados da geometria
        # if plan_type == 1:
        #     ST = (cr + ct)*b/2
        #     mac = (cr + ct)/2
        # else:
        #     S1 = (cr + cm)*b1/2; S2 = (cm + ct)*b2/2
        #     ST = S1 + S2
        #     mac = (S1/b1 + S2/b2)/2
        
        # Escrever o arquivo de entrada
        name = 'wing_geom_DegenGeom'
        fid = open(name + '.vspaero',"w")  # abrir o arquivo
        if plan_type == 0: # Asa trapezoidal simples
            fid.write('Sref = ' + str(ST) + '\n')
            fid.write('Cref = ' + str(mac) + '\n')
            fid.write('Bref = ' + str(b) + '\n')
        else:
            fid.write('Sref = ' + str(ST) + '\n')
            fid.write('Cref = ' + str(mac) + '\n')
            fid.write('Bref = ' + str(b1 + b2) + '\n')
        fid.write('X_cg = 0 \nY_cg = 0\nZ_cg = 0\n')
        fid.write('Mach = 0 \n')
        fid.write('AoA = ' + str(aoa))
        # for i in range(len(aoa)):
        #     fid.write(str(aoa[i]))
        #     if i != len(aoa)-1:
        #         fid.write(', ')
        fid.write('\n')
        fid.write('Beta = 0\n')
        fid.write('Vinf = ' + str(Vinf) + '\n')
        fid.write('Rho = ' + str(rho) + '\n')
        fid.write('ReCref = ' + str(Re_ref) + '\n')
        fid.write('ClMax = -1\n')
        fid.write('MaxTurningAngle = -1\n')
        fid.write('Symmetry = NO\n')
        fid.write('FarDist = -1\n')
        fid.write('NumWakeNodes = 0\n')
        fid.write('WakeIters = ' + str(wake) + '\n')
        fid.write('NumberOfControlGroups = 0\n')
        fid.write('Preconditioner = Matrix\n')
        fid.write('Karman-Tsien Correction = ' + KTC + '\n')
        fid.write('Stability Type = 0\n')
        fid.close() # fechar arquivo
        
        # Executar a simulação
        os.system('vspaero -omp 4 ' + name) 
        # count = 1
        # while 1:
            # p = multiprocessing.Process(target=execute) # Iniciar a execução do VSPAERO como um processo
            # p.start()
            # print('iniciado')
            # # p.join(wait) # Esperar um tempo
            # sleep(wait)
            # print('tempo esgotado')
            
            # if p.is_alive(): # Se ainda estiver funcionando, fechar e continuar o algoritmo
                # print('Terminar e tentar de novo')
                # # Terminate - may not work if process is stuck for good
                # # p.terminate()
                # # OR Kill - will work for sure, no chance for process to finish nicely however
                # p.kill()
    
                # p.join()
            
            # else:
                # p.join(); break
            
            # count += 1
            # if count > tries: 
                # aero = 'n'
                # print('Número máximo de tentativas excedido - prosseguindo para o próximo indivíduo')
                # return aero
        
        
        # Apagar arquivos irrelevantes
        if os.path.exists(name + '.adb'):
            os.remove(name + '.adb')
        if os.path.exists(name + '.adb.cases'):
            os.remove(name + '.adb.cases')
        if os.path.exists(name + '.fem'):
            os.remove(name + '.fem')
        if os.path.exists(name + '.group.1'):
            os.remove(name + '.group.1')
        if os.path.exists(name + '.history'):
            os.remove(name + '.history')
        if os.path.exists(name + '.lod'):
            os.remove(name + '.lod')
        # if os.path.exists(name + '.polar'):
        #     os.remove(name + '.polar')
        # if os.path.exists(name + '.tkey'):
        #     os.remove(name + '.tkey')
        # if os.path.exists(name + '.tri'):
        #     os.remove(name + '.tri')
        # if os.path.exists(name + '.vspaero'):
        #     os.remove(name + '.vspaero')
        
        # Ler o arquivo de polar
        # results = pd.read_csv(name + '.polar',sep = '\t',header = 0) # o campo header indica 'qual linha será usada como cabeçalho
        results = np.genfromtxt(name + '.polar')
        results = results[1:results.shape[0]][:] # tirar a linha com nans
        
        # Obter e mostrar os coeficientes relevantes
        # aero = np.zeros((len(aoa),4))
        
        aero[P][:] = [results[0][4],results[0][7],results[0][9],results[0][15]]
        return aero
        
        # for i in range(len(aoa)):
        #     print('Condição de voo ' + str(i+1) + ': CL = ' + str(aero[i][0]) + ', CD = ' + str(aero[i][1]) + ', L/D = ' + str(aero[i][2]) + ', CM = ' + str(aero[i][3]) + '\n')
        
        
        
        # # Ler o arquivo .tri e fazer um gráfico da asa ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # # Ler todas as coordenadas
        # coo = []
        # with open(name + '.tri') as read_obj:
        #     csv_reader = csv.reader(read_obj,delimiter=' ')
        #     for row in csv_reader:
        #         coo.append(row)
        
        # # Pegar informações
        # n_node = int(coo[0][0]) # Número de nós
        # n_elem = int(coo[0][1]) # Número de elementos
        
        # # Pegar a matriz de elementos
        # elem = coo[n_node+1:n_node+n_elem+1]
        # for i in range(len(elem)):
        #     for j in range(len(elem[i])):
        #         elem[i][j] = int(elem[i][j]) # Converter todos os números pra inteiros
        
        # # Pegar as coordenadas dos nós
        # node = coo[1:n_node+1]
        # for i in range(len(node)):
        #     for j in range(len(node[i])-1,-1,-1):
        #         if node[i][j] == '': 
        #             del(node[i][j]) # Tirar todos os espaços vazios
        #         else:
        #             node[i][j] = float(node[i][j]) # Converter todos os números pra floats
        
        # # Converter matrizes pro formato numpy
        # elem = np.array(elem)
        # node = np.array(node)
        
        # # Ler apenas as coordenadas dos nós (por algum motivo, é isso o que faz a função genfromtxt)
        # # Se for preciso tomar o resto das informações, usar o código acima 
        # # node = np.genfromtxt(name + '.tri')
        
        # # Fazer o gráfico
        # fig1 = plt.figure()
        # ax1 = fig1.add_subplot(projection='3d')
        # ax1.scatter(node[:,0],node[:,1],node[:,2]) # Mostrar os nós
        # ax1.set_xlabel('X Label')
        # ax1.set_ylabel('Y Label')
        # ax1.set_zlabel('Z Label')
        # # Create cubic bounding box to simulate equal aspect ratio
        # max_range = np.array([node[:,0].max()-node[:,0].min(), node[:,1].max()-node[:,1].min(), node[:,2].max()-node[:,2].min()]).max()
        # Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(node[:,0].max()+node[:,0].min())
        # Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(node[:,1].max()+node[:,1].min())
        # Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(node[:,2].max()+node[:,2].min())
        # # Comment or uncomment following both lines to test the fake bounding box:
        # for xb, yb, zb in zip(Xb, Yb, Zb):
        #    ax1.plot([xb], [yb], [zb], 'w')
        # ax1.view_init(90,90) # usar isto pra gerar um pop-up da janela do gráfico
        
        # # Montar os elementos
        # for i in range(len(elem)):
        #     n1 = elem[i,0]
        #     n2 = elem[i,1]
        #     n3 = elem[i,2]
        #     ax1.plot([node[n1,0],node[n1,1],node[n1,2]],[node[n2,0],node[n2,1],node[n2,2]])
        #     ax1.plot([node[n2,0],node[n2,1],node[n2,2]],[node[n3,0],node[n3,1],node[n3,2]])
        #     ax1.plot([node[n3,0],node[n3,1],node[n3,2]],[node[n1,0],node[n1,1],node[n1,2]])
        
        
        # p1 = np.array([1,1])
        # p2 = np.array([2,2])
        # plt.plot(,)
        
        
        
