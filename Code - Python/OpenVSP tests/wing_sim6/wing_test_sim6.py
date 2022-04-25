# Objetivo: gerar uma asa e simular
# Descrição desta versão: o código gera uma asa trapezoidal simples ou dupla,
# de acordo com a configuração. É possível alterar a planta da asa e a torção.
# Os aerofólios são naca 4 dígitos. Ao final é feita uma simulação no prompt de
# comando, os resultados são lidos, e é apresentada uma figura da asa.

# A simulação é feita pelo VLM, pela API (ao invés de invocar o CMD)


 
import openvsp as vsp
import math
import os
import numpy as np
import csv
import matplotlib.pyplot as plt


name = 'wing_test_sim6'

# Reiniciar o modelo 
vsp.ClearVSPModel()

# // ===== Criar uma asa trapezoidal simples ou dupla com perfis NACA 4 dígitos ====//
# Nota: todas as unidades no sistema internacional

# Dados da geometria
plan_op = 1 # Opção de forma da planta -> 1 para trapezoidal simples, outro valor para trapezoidal dupla
# ST = 30 # Área alar (total - lembrar que o openvsp trata os valores para cada seção independentemente)
# S1 = 15 # Área alar da primeira seção da asa (caso a asa seja trapezoidal dupla)
# S2 = ST - S1 # Área da segunda seção da asa (caso a asa seja trapezoidal dupla)
b = 5
b1 = 15
b2 = b - b1
cr = 3; # Corda da raiz
cm = 2; # Corda do meio
ct = 1; # Corda da ponta
naca_r = '2415' # perfil da raiz
naca_m = '2412' # perfil do meio
naca_t = '0012' # perfil da ponta
tw_m = -0 # Torção geométrica no meio
tw_t = -0 # Torção geométrica na ponta
# Parâmetros de qualidade da malha (mosaico) - isto afeta diretamente o arquivo .tri, o que por sua vez afeta os resultados das simulações
num_U1 = 20
num_W1 = 40
num_U1 = 20
num_W1 = 40


# Dados da simulação
sim_op = 0 # Atribuir valor 1 se quiser fazer a simulação
aoa = [0,2,4] # Ângulos de ataque
Vinf = 100 # Velocidade
Re_ref = 1e6 # Número de Reynolds na corda de referência


# Montar a geometria
if plan_op == 1: # Asa trapezoidal simples
    # Calcular enflechamento do bordo de ataque, considerando que a linha c/2 não
    # deve ter enflechamento
    # b = ST*2/(cr + ct)
    lambda_le = math.degrees(math.atan((cr - ct)/b))

    #//==== Add Wing ====//
    wid = vsp.AddGeom( "WING", "" )

    #//===== Change Driver =====//
    # isto funciona da mesma forma que no gui: o projetista determina três parâmetros
    # entre os 8 possíveis pra manipular a planta da asa
    # vsp.SetDriverGroup( wid, 1, vsp.AREA_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    vsp.SetDriverGroup( wid, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )

    # vsp.SetParmVal( wid, "RotateAirfoilMatchDideralFlag", "WingGeom", 1.0 )

    #//===== Change Some Parameters 1st Section ====//
    vsp.SetParmVal( wid, "Root_Chord", "XSec_1", cr )
    vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", ct )
    # vsp.SetParmVal( wid, "Area", "XSec_1", ST/2 )
    vsp.SetParmVal( wid, "Span", "XSec_1", b/2 )
    vsp.SetParmVal( wid, "Sweep", "XSec_1",lambda_le )
    if tw_t != 0: vsp.SetParmVal(wid,"Twist","XSec_1",tw_t)
    vsp.Update()

    # Mudar aerofólio da raiz
    vsp.SetParmVal( wid, "Camber", "XSecCurve_0", float(naca_r[0])/100 ) # Curvatura máxima (m)
    vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_0", float(naca_r[1])/10 )  # Local da curvatura máxima (p)
    vsp.SetParmVal( wid, "ThickChord", "XSecCurve_0", float(naca_r[2:4])/100 )  # Espessura máxima (t)
    vsp.Update()

    # Mudar aerofólio da ponta
    vsp.SetParmVal( wid, "Camber", "XSecCurve_1", float(naca_t[0])/100 ) # Curvatura máxima (m)
    vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_1", float(naca_t[1])/10 )  # Local da curvatura máxima (p)
    vsp.SetParmVal( wid, "ThickChord", "XSecCurve_1", float(naca_t[2:4])/100 )  # Espessura máxima (t)
    vsp.Update()
    
    # Alterar parâmetros de qualidade
    if num_W1 != 0: vsp.SetParmVal(wid,'Tess_W','Shape',num_W1)
    if num_U1 != 0: vsp.SetParmVal(wid,'SectTess_U','XSec_1',num_U1)
    

else: # Asa trapezoidal dupla
    # Calcular enflechamento do bordo de ataque da primeira seção
    # b1 = S1*2/(cr + cm)
    # b2 = S2*2/(cm + ct)
    lambda_le1 = math.degrees(math.atan((cr - cm)/b1))
    lambda_le2 = math.degrees(math.atan((cm - ct)/b2))

    #//==== Add Wing ====//
    wid = vsp.AddGeom( "WING", "" )
    
    # Adicionar mais uma seção
    vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )

    #//===== Change Driver =====//
    # isto funciona da mesma forma que no gui: o projetista determina três parâmetros
    # entre os 8 possíveis pra manipular a planta da asa
    vsp.SetDriverGroup( wid, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    vsp.SetDriverGroup( wid, 2, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )

    # vsp.SetParmVal( wid, "RotateAirfoilMatchDideralFlag", "WingGeom", 1.0 )

    #//===== Change Some Parameters 1st Section ====//
    vsp.SetParmVal( wid, "Root_Chord", "XSec_1", cr )
    vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", cm )
    # vsp.SetParmVal( wid, "Area", "XSec_1", S1/2 )
    vsp.SetParmVal( wid, "Span", "XSec_1", b1/2 )
    vsp.SetParmVal( wid, "Sweep", "XSec_1",lambda_le1 )
    if tw_t != 0: vsp.SetParmVal(wid,"Twist","XSec_1",tw_m)
    vsp.Update()
    
    # Mudar parâmetros da segunda seção da asa
    # vsp.SetParmVal( wid, "Root_Chord", "XSec_2", cr )
    vsp.SetParmVal( wid, "Tip_Chord", "XSec_2", ct )
    # vsp.SetParmVal( wid, "Area", "XSec_2", S2/2 )
    vsp.SetParmVal( wid, "Span", "XSec_2", b2/2 )
    vsp.SetParmVal( wid, "Sweep", "XSec_2",lambda_le2 )
    if tw_t != 0: vsp.SetParmVal(wid,"Twist","XSec_2",tw_t)
    vsp.Update()

    # Mudar aerofólio da raiz
    vsp.SetParmVal( wid, "Camber", "XSecCurve_0", float(naca_r[0])/100 ) # Curvatura máxima (m)
    vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_0", float(naca_r[1])/10 )  # Local da curvatura máxima (p)
    vsp.SetParmVal( wid, "ThickChord", "XSecCurve_0", float(naca_r[2:4])/100 )  # Espessura máxima (t)
    vsp.Update()

    # Mudar aerofólio do meio
    vsp.SetParmVal( wid, "Camber", "XSecCurve_1", float(naca_m[0])/100 ) # Curvatura máxima (m)
    vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_1", float(naca_m[1])/10 )  # Local da curvatura máxima (p)
    vsp.SetParmVal( wid, "ThickChord", "XSecCurve_1", float(naca_m[2:4])/100 )  # Espessura máxima (t)
    vsp.Update()

    # Mudar aerofólio da ponta
    vsp.SetParmVal( wid, "Camber", "XSecCurve_2", float(naca_t[0])/100 ) # Curvatura máxima (m)
    vsp.SetParmVal( wid, "CamberLoc", "XSecCurve_2", float(naca_t[1])/10 )  # Local da curvatura máxima (p)
    vsp.SetParmVal( wid, "ThickChord", "XSecCurve_2", float(naca_t[2:4])/100 )  # Espessura máxima (t)
    vsp.Update()
    
    
    
# //==== Check For API Errors ====//
# while ( GetNumTotalErrors() > 0 )
# {
#     ErrorObj err = PopLastError();
#     Print( err.GetErrorString() );
#    }
errorMgr = vsp.ErrorMgrSingleton_getInstance()
num_err = errorMgr.GetNumTotalErrors()
for i in range(0, num_err):
    err = errorMgr.PopLastError()
    print("error = ", err.m_ErrorString)

# Exportar arquivo do openvsp
vsp.WriteVSPFile(name + ".vsp3")

# Salvar a malha .tri do modelo
# vsp.ExportFile(name + '.tri',vsp.SET_ALL,vsp.EXPORT_CART3D )



# Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if plan_op == 1:
    ST = (cr + ct)*b/2
    mac = (cr + ct)/2
else:
    S1 = (cr + cm)*b1/2; S2 = (cm + ct)*b2/2
    ST = S1 + S2
    mac = (S1/b1 + S2/b2)/2

if sim_op == 1:
    # Escrever o arquivo de entrada
    fid = open(name + '.vspaero',"w")  # abrir o arquivo
    if plan_op == 1: # Asa trapezoidal simples
        fid.write('Sref = ' + str(ST) + '\n')
        fid.write('Cref = ' + str((cr + ct)/2) + '\n')
        fid.write('Bref = ' + str(b) + '\n')
    else:
        fid.write('Sref = ' + str(ST) + '\n')
        fid.write('Cref = ' + str((cr + cm + ct)/3) + '\n')
        fid.write('Bref = ' + str(b1 + b2) + '\n')
    fid.write('X_cg = 0 \nY_cg = 0\nZ_cg = 0\n')
    fid.write('Mach = 0 \n')
    fid.write('AoA = ')
    for i in range(len(aoa)):
        fid.write(str(aoa[i]))
        if i != len(aoa)-1:
            fid.write(', ')
    fid.write('\n')
    fid.write('Beta = 0\n')
    fid.write('Vinf = ' + str(Vinf) + '\n')
    fid.write('Rho = 0.002377\n')
    fid.write('ReCref = 10000000\n')
    fid.write('ClMax = -1\n')
    fid.write('MaxTurningAngle = -1\n')
    fid.write('Symmetry = NO\n')
    fid.write('FarDist = -1\n')
    fid.write('NumWakeNodes = 0\n')
    fid.write('WakeIters = 5\n')
    fid.write('NumberOfControlGroups = 0\n')
    fid.write('Preconditioner = Matrix\n')
    fid.write('Karman-Tsien Correction = N\n')
    fid.write('Stability Type = 0\n')
    fid.close() # fechar arquivo
    
    # Executar a simulação
    os.system('vspaero ' + name) 
    
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
    if os.path.exists(name + '.tkey'):
        os.remove(name + '.tkey')
    # if os.path.exists(name + '.tri'):
    #     os.remove(name + '.tri')
    # if os.path.exists(name + '.vspaero'):
    #     os.remove(name + '.vspaero')
    
#     # Ler o arquivo de polar
#     # results = pd.read_csv(name + '.polar',sep = '\t',header = 0) # o campo header indica 'qual linha será usada como cabeçalho
#     results = np.genfromtxt(name + '.polar')
#     results = results[1:results.shape[0]][:] # tirar a linha com nans
    
#     # Obter e mostrar os coeficientes relevantes
#     aero = np.zeros((len(aoa),4))
#     for i in range(len(aoa)):
#         aero[i][:] = [results[i][4],results[i][7],results[i][9],results[i][15]]
    
#     for i in range(len(aoa)):
#         print('Condição de voo ' + str(i+1) + ': CL = ' + str(aero[i][0]) + ', CD = ' + str(aero[i][1]) + ', L/D = ' + str(aero[i][2]) + ', CM = ' + str(aero[i][3]) + '\n')



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
# # ax1.view_init(90,90) # usar isto pra gerar um pop-up da janela do gráfico

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





# vsp.ComputeDegenGeom( vsp.SET_ALL, vsp.DEGEN_GEOM_CSV_TYPE )