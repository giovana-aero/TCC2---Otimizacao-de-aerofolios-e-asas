# Objetivo: gerar a geometria de uma asa no openvsp, cujos aerofólios
# são definidos por arquivos de coordenadas

import openvsp as vsp
import numpy as np
import math
import csv
import matplotlib.pyplot as plt




# Reiniciar o modelo 
vsp.ClearVSPModel()


# // ===== Criar uma asa trapezoidal simples ====//
# Nota: todas as unidades no sistema internacional

# Dados da geometria
S = 30; # Área alar (total - lembrar que o openvsp trata os valores para cada seção independentemente)
cr = 3; # Corda da raiz
ct = 1; # Corda da ponta
af_r = 'Wortmann FX 61-184.dat' # perfil da raiz
af_t = 'Wortmann FX 61-184.dat' # perfil da ponta

# Calcular enflechamento do bordo de ataque, considerando que a linha c/2 não
# deve ter enflechamento
b = S*2/(cr + ct)
lambda_le = math.degrees(math.asin((cr - ct)/b))


#//==== Add Wing ====//
wid = vsp.AddGeom( "WING", "" )


#//===== Change Driver =====//
# isto funciona da mesma forma que no gui: o projetista determina três parâmetros
# entre os 8 possíveis pra manipular a planta da asa
vsp.SetDriverGroup( wid, 1, vsp.AREA_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )

# vsp.SetParmVal( wid, "RotateAirfoilMatchDideralFlag", "WingGeom", 1.0 )

#//===== Change Some Parameters 1st Section ====//
vsp.SetParmVal( wid, "Root_Chord", "XSec_1", cr )
vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", ct )
vsp.SetParmVal( wid, "Area", "XSec_1", S/2 )
vsp.SetParmVal( wid, "Sweep", "XSec_1",lambda_le )
vsp.Update()

# Mudar aerofólio da raiz
xsec_surf = vsp.GetXSecSurf( wid, 0 ) # Isto vale pra asa completa
vsp.ChangeXSecShape( xsec_surf, 0, vsp.XS_FILE_AIRFOIL )
xsec_r = vsp.GetXSec( xsec_surf, 0 )
vsp.ReadFileAirfoil( xsec_r, af_r )
vsp.Update()

# Mudar aerofólio da ponta
vsp.ChangeXSecShape( xsec_surf, 1, vsp.XS_FILE_AIRFOIL )
xsec_t = vsp.GetXSec( xsec_surf, 1 )
vsp.ReadFileAirfoil( xsec_t, af_t )
vsp.Update()

# Salvar modelo
vsp.WriteVSPFile("wing_airfoil_coo.vsp3")

# Salvar malha
vsp.ExportFile('wing_airfoil_coo.tri',vsp.SET_ALL,vsp.EXPORT_CART3D )

# Ler o arquivo .tri e fazer um gráfico da asa ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ler todas as coordenadas
coo = []
with open('wing_airfoil_coo.tri') as read_obj:
    csv_reader = csv.reader(read_obj,delimiter=' ')
    for row in csv_reader:
        coo.append(row)

# Pegar informações
n_node = int(coo[0][0]) # Número de nós
n_elem = int(coo[0][1]) # Número de elementos

# Pegar a matriz de elementos
elem = coo[n_node+1:n_node+n_elem+1]
for i in range(len(elem)):
    for j in range(len(elem[i])):
        elem[i][j] = int(elem[i][j]) # Converter todos os números pra inteiros

# Pegar as coordenadas dos nós
node = coo[1:n_node+1]
for i in range(len(node)):
    for j in range(len(node[i])-1,-1,-1):
        if node[i][j] == '': 
            del(node[i][j]) # Tirar todos os espaços vazios
        else:
            node[i][j] = float(node[i][j]) # Converter todos os números pra floats

# Converter matrizes pro formato numpy
elem = np.array(elem)
node = np.array(node)

# Fazer o gráfico
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection='3d')
ax1.scatter(node[:,0],node[:,1],node[:,2]) # Mostrar os nós
ax1.set_xlabel('X Label')
ax1.set_ylabel('Y Label')
ax1.set_zlabel('Z Label')
# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array([node[:,0].max()-node[:,0].min(), node[:,1].max()-node[:,1].min(), node[:,2].max()-node[:,2].min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(node[:,0].max()+node[:,0].min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(node[:,1].max()+node[:,1].min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(node[:,2].max()+node[:,2].min())
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax1.plot([xb], [yb], [zb], 'w')
ax1.view_init(45,45) # usar isto pra gerar um pop-up da janela do gráfico (na verdade isto não funciona pra isso?)




## Anotações

# NOTA: o carregamento de arquivos de aerofólio pode ter algo a ver com o enum XS_FILE_AIRFOIL 
# ver http://openvsp.org/api_docs/latest/group___enumerations.html

# Descrições das entradas das funções
# xsec_surf = vsp.GetXsecSurf(geom_id,XSecSurf index)
# vsp.ChangeXSecShape(XsecSurf id, Xsec index, Xsec type enum)
# xsec = vsp.GetXSec(XsecSurf id, Xsec index)
# vso.ReadFileAirfoil(Xsec id, file_name)



# Teste de fuselagem (visto na documentação da função ReadFileAirfoil)
# vsp.ClearVSPModel()
# fuseid =vsp.AddGeom( "FUSELAGE", "" );
# xsec_surf = vsp.GetXSecSurf( fuseid, 0 );
# vsp.ChangeXSecShape( xsec_surf, 1, vsp.XS_FILE_AIRFOIL );
# xsec = vsp.GetXSec( xsec_surf, 1 ); 
# vsp.ReadFileAirfoil( xsec, "coordenadas.dat" );
# vsp.Update()
# vsp.WriteVSPFile( "aoba.vsp3")


# Teste da asa, baseado no teste da fuselagem
# vsp.ClearVSPModel()
# # wid = vsp.AddGeom('WING')
# xsec_surf = vsp.GetXSecSurf( wid, 0 ) # Isto vale pra asa completa
# vsp.ChangeXSecShape( xsec_surf, 0, vsp.XS_FILE_AIRFOIL )
# xsec_r = vsp.GetXSec( xsec_surf, 0 )
# vsp.ReadFileAirfoil( xsec_r, af_r )
# vsp.Update()
# vsp.WriteVSPFile( "aoba.vsp3")
