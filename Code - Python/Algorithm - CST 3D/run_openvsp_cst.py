import openvsp as vsp
import math
import os
from airfoil_interpolation import airfoil_interpolation
from run_cst_TCC2_3D import run_cst_TCC2_3D
import numpy as np


def run_openvsp_cst(pop,dat):
    
    # Criar a geometria da asa e escrever o arquivo .vsp3 ~~~~~~~~~~~~~~~~~~~~~
    
    name = 'wing_geom'
    
    # Reiniciar o modelo 
    vsp.ClearVSPModel()
    
    # Dados da geometria
    plan_type = pop.type # Opção de forma da planta -> 0 para trapezoidal simples, outro valor para trapezoidal dupla
    b = pop.b # Envergadura total
    if plan_type == 1: 
        b1 = pop.b1 # Envergadura da primeira seção
        b2 = b - b1 # Envergadura da segunda seção
    cr = pop.c_r # Corda da raiz
    cm = pop.c_m # Corda do meio
    ct = pop.c_t # Corda da ponta
    if plan_type == 1: tw_m = pop.tw_m # Torção geométrica no meio
    tw_t = pop.tw_t # Torção geométrica na ponta
    # Parâmetros de qualidade da malha (mosaico) - isto afeta diretamente o arquivo .tri, o que por sua vez afeta os resultados das simulações
    # num_U1 = dat.num_U1
    # num_W1 = dat.num_W1
    
    # Gerar as coordenadas dos perfis da raiz e da ponta
    run_cst_TCC2_3D(pop.v_ex_r,pop.v_in_r,dat,[dat.N1_r,dat.N2_r],1,'coo_r.dat') # Raiz
    run_cst_TCC2_3D(pop.v_ex_t,pop.v_in_t,dat,[dat.N1_t,dat.N2_t],1,'coo_t.dat') # Ponta
    
    # Montar a geometria
    if plan_type == 0: # Asa trapezoidal simples
        # Calcular enflechamento do bordo de ataque, considerando que a linha c/2 não
        # deve ter enflechamento
        # b = ST*2/(cr + ct)
        lambda_le = math.degrees(math.atan((cr - ct)/b))
    
        # Adicionar asa
        wid = vsp.AddGeom( "WING", "" )
    
        # Mudar drivers
        # isto funciona da mesma forma que no gui: o projetista determina três parâmetros
        # entre os 8 possíveis pra manipular a planta da asa
        vsp.SetDriverGroup( wid, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    
        # Mudar parâmetros da planta
        vsp.SetParmVal( wid, "Root_Chord", "XSec_1", cr ) # Corda da raiz
        vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", ct ) # Corda da ponta
        vsp.SetParmVal( wid, "Span", "XSec_1", b/2 ) # Envergadura
        vsp.SetParmVal( wid, "Sweep", "XSec_1",lambda_le ) # Enflechamento
        if tw_t != 0: vsp.SetParmVal(wid,"Twist","XSec_1",tw_t) # Torção geométrica
        vsp.Update() 
    
        # Mudar aerofólio da raiz
        xsec_surf = vsp.GetXSecSurf( wid, 0 ) # Isto vale pra asa completa
        vsp.ChangeXSecShape( xsec_surf, 0, vsp.XS_FILE_AIRFOIL )
        xsec_r = vsp.GetXSec( xsec_surf, 0 )
        vsp.ReadFileAirfoil( xsec_r, 'coo_r.dat' )
        vsp.Update()
    
        # Mudar aerofólio da ponta
        vsp.ChangeXSecShape( xsec_surf, 1, vsp.XS_FILE_AIRFOIL )
        xsec_t = vsp.GetXSec( xsec_surf, 1 )
        vsp.ReadFileAirfoil( xsec_t, 'coo_t.dat' )
        vsp.Update()
        
        # Alterar parâmetros de qualidade da malha
        # if num_W1 != 0: vsp.SetParmVal(wid,'Tess_W','Shape',num_W1)
        if dat.nb[0] != 0: 
            # Alterar a quantidade de seções em termos de uma concentração especificada
            if dat.nb[1] == 0: # Usar o valor especificado 
                nb = dat.nb[0]
            else: # Usar como uma concentração por metro e determinar o número de seções
                nb = np.floor((b*dat.nb[0]-2)/2) + 2
            
            vsp.SetParmVal(wid,'SectTess_U','XSec_1',nb)
        
    
    else: # Asa trapezoidal dupla
        
        # Considerar a opção 'L' da corda do meio
        if cm == 'L':
            cm = 2*(ct-cr)/b*b1/2 + cr
            
        # Considerar a opção 'L' do aerofólio do meio
        if 'L' in dat.le_R_ext1_in_m:
            coo1 = np.genfromtxt('coo_r.dat',delimiter=' ')
            coo2 = np.genfromtxt('coo_t.dat',delimiter=' ')
            coo_intp = airfoil_interpolation(coo1,coo2,b1/b)
            with open('coo_m.dat','w')as f: # Imprimir pra ser lido posteriormente 
                np.savetxt(f, coo_intp, delimiter=' ', newline='\n', header='', footer='', comments='# ')
        else: # Gerar as coordenadas normalmente
            run_cst_TCC2_3D(pop.v_ex_m,pop.v_in_m,dat,[dat.N1_r,dat.N2_r],1,'coo_m.dat')
                
        # Considerar a opção 'L' da torção geométrica do meio e da ponta
        if tw_m == 'L' and not isinstance(tw_t,str): # Definir a torção do meio em termos da torção da ponta
            tw_m = tw_t*b1/b
        elif tw_t == 'L' and not isinstance(tw_m,str): # Definir a torção da ponta em termos da torção do meio
            tw_t = tw_m*b/b1

        
        
    
    
        # Calcular enflechamentos dos bordos de ataque 
        lambda_le1 = math.degrees(math.atan((cr - cm)/b1))
        lambda_le2 = math.degrees(math.atan((cm - ct)/b2))
    
        # Adicionar asa
        wid = vsp.AddGeom( "WING", "" )
        
        # Adicionar mais uma seção
        vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )
    
        # Mudar drivers
        # isto funciona da mesma forma que no gui: o projetista determina três parâmetros
        # entre os 8 possíveis pra manipular a planta da asa
        vsp.SetDriverGroup( wid, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
        vsp.SetDriverGroup( wid, 2, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    
        # Mudar parâmetros da planta (primeira seção)
        vsp.SetParmVal( wid, "Root_Chord", "XSec_1", cr ) # Corda da raiz
        vsp.SetParmVal( wid, "Tip_Chord", "XSec_1", cm ) # Corda do meio (ponta desta seção)
        vsp.SetParmVal( wid, "Span", "XSec_1", b1/2 ) # Envergadura
        vsp.SetParmVal( wid, "Sweep", "XSec_1",lambda_le1 ) # Enflechamento
        if tw_m != 0: vsp.SetParmVal(wid,"Twist","XSec_1",tw_m) # Torção geométrica
        vsp.Update()
        
        # Mudar parâmetros da planta (segunda seção)
        vsp.SetParmVal( wid, "Tip_Chord", "XSec_2", ct ) # Corda da ponta
        vsp.SetParmVal( wid, "Span", "XSec_2", b2/2 ) # Envergadura
        vsp.SetParmVal( wid, "Sweep", "XSec_2",lambda_le2 ) # Enflechamento
        if tw_t != 0: vsp.SetParmVal(wid,"Twist","XSec_2",tw_t) # Torção geométrica
        vsp.Update()
    
        # Mudar aerofólio da raiz
        xsec_surf = vsp.GetXSecSurf( wid, 0 ) # Isto vale pra asa completa
        vsp.ChangeXSecShape( xsec_surf, 0, vsp.XS_FILE_AIRFOIL )
        xsec_r = vsp.GetXSec( xsec_surf, 0 )
        vsp.ReadFileAirfoil( xsec_r, 'coo_r.dat' )
        vsp.Update()

        # Mudar aerofólio do meio
        # if 'L' in dat.le_R_ext1_in_m: # Considerar a opção 'L' do aerofólio do meio
        #     xsec_surf = vsp.GetXSecSurf( wid, 0 )     
        #     vsp.ChangeXSecShape( xsec_surf, 1, vsp.XS_FILE_AIRFOIL )
        #     xsec_m = vsp.GetXSec( xsec_surf, 1 )
        #     vsp.ReadFileAirfoil( xsec_m, 'coo_m.dat'  )
        #     vsp.Update()    

        # else: # Inserir as informações convencionalmente
        vsp.ChangeXSecShape( xsec_surf, 1, vsp.XS_FILE_AIRFOIL )
        xsec_m = vsp.GetXSec( xsec_surf, 1 )
        vsp.ReadFileAirfoil( xsec_m, 'coo_m.dat' )
        vsp.Update()
    
        # Mudar aerofólio da ponta
        vsp.ChangeXSecShape( xsec_surf, 2, vsp.XS_FILE_AIRFOIL )
        xsec_t = vsp.GetXSec( xsec_surf, 2 )
        vsp.ReadFileAirfoil( xsec_t, 'coo_t.dat' )
        vsp.Update()
    
        # Alterar parâmetros de qualidade da malha
        # if num_W1 != 0: vsp.SetParmVal(wid,'Tess_W','Shape',num_W1)
        if dat.nb[0] != 0: 
            if isinstance(dat.nb1,str) and dat.nb1 == 'L': # Tornar uniforme (ou próximo disso) a distribuição de seções ao longo da envergadura completa
                nb = dat.nb
                # Alterar a quantidade de seções em termos de uma concentração especificada
                if nb[1] == 0: # Usar o valor especificado 
                    nb = nb[0] + 2
                else: # Usar como uma concentração por metro e determinar o número de seções
                    nb = np.floor((b*nb[0]-2)/2) + 2
                
                nb1 = nb*b1/b
                if nb1 - round(nb1) < 0: # Caso o valor seja arredondado para cima
                    nb1 = round(nb1)
                    nb2 = nb - nb1
                else: # Caso o valor seja arredondado para baixo
                    nb1 = round(nb1)
                    nb2 = nb - nb1
                
            else: # Usar os valores originais de nb1 e nb2 dados pelo usuário
                nb2 = dat.nb2
                nb = nb1 + nb2
            
            vsp.SetParmVal(wid,'SectTess_U','XSec_1',nb1)
            vsp.SetParmVal(wid,'SectTess_U','XSec_2',nb2)
        
        
    # Checar erros na API
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
    vsp.ExportFile(name + '.tri',vsp.SET_ALL,vsp.EXPORT_CART3D )
    os.remove(name + '.tkey')
    
    # Mais alguns dados da geometria
    if plan_type == 0:
        ST = (cr + ct)*b/2
        mac = (cr + ct)/2
    else:
        S1 = (cr + cm)*b1/2; S2 = (cm + ct)*b2/2
        ST = S1 + S2
        mac = (S1/b1 + S2/b2)/2
    pop.mac = mac
    pop.S = ST
    
    return pop
    
    