import numpy as np
from math import radians,degrees,cos,tan,atan
from ypstruct import structure
from airfoil_interpolation import airfoil_interpolation
from airfoil_rotation import airfoil_rotation
from numpy.matlib import repmat
from run_cst_TCC2_3D import run_cst_TCC2_3D
import matplotlib.pyplot as plt

def run_apame_mesher_cst(pop,dat,op=0):
    # Esta função gera a malha de asas pro apame. Asas são trapezoidais simples ou
    # duplas.
    
    # Nota sobre as coordenadas de aerofólios: este script trabalha de modo que seja
    # necessário que os aerofólios carregados tenham sempre o mesmo número de pontos.
    # Isso não é um problema no contexto do algoritmo de otimização. No entanto, se o
    # o intuito for gerar uma malha em outra aplicação, é necessário garantir esse
    # requisito de números iguais de pontos. Para isso, pode-se interpolar as coordenadas
    # em mãos, ou inserí-las no CST reverso e gerar o aerofólio com a quantidade desejada
    # de pontos
    
    # Valores de entrada para asas trapezoidais simples:
    # b     - envergadura
    # c_r   - corda da raiz
    # c_t   - corda da ponta
    # tw_t  - torção geométrica na ponta
    # af_r - nome do aerofólio da raiz [m,p,t]
    # af_t - nome do aerofólio da ponta [m,p,t]
    # nb    - número de seções intermediárias (raiz/ponta)
    
    # Valores de entrada para asas trapezoidais duplas:
    # b    - envergadura
    # b1   - envergadura da primeira parte da asa (raiz ao meio)
    # c_r  - corda da raiz
    # c_m  - corda do meio (permite opção 'L')
    # c_t  - corda da ponta
    # tw_m - torção geométrica no meio (permite opção 'L')
    # tw_t - torção geométrica na ponta
    # af_r - nome do aerofólio da raiz [m,p,t]
    # af_m - nome do aerofólio do meio [m,p,t] (permite opção 'L')
    # af_t - nome do aerofólio da ponta [m,p,t]
    # nb   - número de seções intermediárias (raiz/ponta) 
    # nb1  - número de seçoes intermediárias (raiz/meio) (permite opção 'L')
    # nb2  - número de seçoes intermediárias (meio/ponta)
    
    # Para visualizar um gráfico da asa:
    # op1 = 2 -> vista da planta
    # op1 = 3 -> vista isométrica
    # Para visualizar os aerofólios principais: 
    # op2 = 1 -> formatos originais
    # op2 = 2 -> formatos com rotação
    
    if pop.type == 0: # Se a asa for trapezoidal simples

        # Dados da asa
        b = pop.b # Envergadura
        c_r = pop.c_r # Corda da raiz
        c_t = pop.c_t # Corda da ponta (DEVE ser menor que a da raiz)
        tw_t = pop.tw_t # Torção geométrica na ponta
        sweep = pop.sweep # Enflechamento total

        # Carregar coordenadas dos aerofólios. Como o contorno é fechado, ignora-se o último par de coordenadas
        coo_r = run_cst_TCC2_3D(pop.v_ex_r,pop.v_in_r,dat,[dat.N1_r,dat.N2_r]); coo_r = coo_r[0:-1,:]*c_r
        coo_t = run_cst_TCC2_3D(pop.v_ex_t,pop.v_in_t,dat,[dat.N1_t,dat.N2_t]); coo_t = coo_t[0:-1,:]*c_t

        # Dados da malha de painéis
        nb = dat.nb
        far = b*2 # Comprimento dos painéis de trilha (a partir do bordo de fuga)
        
        # Alterar a quantidade de seções em termos de uma concentração especificada
        if nb[1] == 0: # Usar o valor especificado 
            nb = nb[0]
        else: # Usar como uma concentração por metro e determinar o número de seções
            nb = int(np.floor((b*nb[0]-2)/2))
        
        # Alguns dados a mais
        sec_af_N = coo_r.shape[0] # Número de nós por seção
        sec_N = 3 + nb*2 # Número de seções transversais

        # Rotacionar o perfil da ponta
        coo_t_R = airfoil_rotation(coo_t,tw_t)
        coo_t[:,0] += 1/4 # Consertar as coordenadas (coo_t é afetado pela função airfoil_rotation. não sei por quê)

        # Obter nós ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        NODE = np.zeros((sec_af_N*(1+nb),3))
        if sweep == 'Z': # Fazer com que o enflechamento da linha c/2 seja sempre zero
            dx_t = (c_r - c_t*cos(radians(tw_t)))/2 # Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
        else: # Usar enflechamento especificado pelo usuário
            dx_t = b/2*tan(radians(sweep))
        
        NODE[0:sec_af_N,:] = np.hstack((np.array([coo_t_R[:,0]]).transpose()+dx_t,repmat(-b/2,sec_af_N,1),np.array([coo_t_R[:,1]]).transpose())) # Ponta esquerda da asa

        # Gerar seções intermediárias e adicionar ao lado esquerdo da asa
        if nb > 0:
            # Criar struct que guarda seções de asa intermediárias
            # (isto será um auxílio devido ànatureza simétrica a asa)
            wing_sec = structure()
            wing_sec.coo = None; wing_sec.dx = 0
            wing_sec = wing_sec.repeat(nb) # Inicializar
            op_vec = np.linspace(1,0,2+nb); op_vec = op_vec[1:-1] # Definição do formato da interpolação em função dos originais
            
            # Encontrar as coordenadas das seções (fazer interpolações)
            for i in range(len(op_vec)):
                # Fazer interpolação e aplicar torção geométrica
                tw_sec = op_vec[i]*tw_t
                wing_sec[i].coo = airfoil_interpolation(coo_r,coo_t,op_vec[i],-op_vec[i]*(b/2),tw_sec)
                # Aplicar o translado
                wing_sec[i].dx = op_vec[i]*dx_t
                wing_sec[i].coo[:,0] += wing_sec[i].dx
                # Adicionar ao lado esquerdo da asa
                NODE[sec_af_N*(i+1):sec_af_N*(i+2),:] = wing_sec[i].coo

        # Adicionar coordenadas do perfil da raiz
        NODE = np.vstack((NODE,   np.hstack((np.array([coo_r[:,0]]).transpose(),np.zeros((sec_af_N,1)),np.array([coo_r[:,1]]).transpose()))))

        # Adicionar seções intermediárias ao lado direito da asa
        if nb > 0:
            temp = np.zeros((sec_af_N*nb,3))
            k = 0
            for i in range(len(op_vec)-1,-1,-1):
                wing_sec[k].coo[:,1] = -wing_sec[k].coo[:,1] # Inverter o sinal da coordenada y das seções intermediárias
                temp[sec_af_N*(i+1)-sec_af_N:sec_af_N*(i+1),:] = wing_sec[k].coo
                k += 1
            
            NODE = np.vstack((NODE,temp))
        
        # Adicionar ponta direita
        NODE = np.vstack((NODE,np.hstack((np.array([coo_t_R[:,0]]).transpose()+dx_t,np.zeros((sec_af_N,1))+b/2,np.array([coo_t_R[:,1]]).transpose()))))

        # Calcular área e corda aerodinâmica média
        S = (c_r + c_t)*b/2
        mac = (c_r + c_t)/2

    else:# Asa trapezoidal dupla
        # Dados da asa
        b = pop.b # Envergadura
        b1 = pop.b1 # Envergadura da primeira seção (raiz ao meio)
        b2 = b - b1 # Envergadura da segunda seção (meio à raiz)
        c_r = pop.c_r # Corda da raiz
        c_m = pop.c_m # Corda do meio
        c_t = pop.c_t # Corda da ponta (DEVE ser menor que a da raiz)
        tw_m = pop.tw_m # Torção geométrica no meio
        tw_t = pop.tw_t # Torção geométrica na ponta
        sweep1 = pop.sweep1 # Enflechamento da primeira seção
        sweep2 = pop.sweep2 # Enflechamento da segunda seção
        
        # Considerar opção 'L' da corda do meio
        if c_m == 'L':
            c_m = 2*(c_t-c_r)/b*b1/2 + c_r
        
        # Carregar coordenadas dos aerofólios da raiz e da ponta. Como o contorno é fechado, ignora-se o último par de coordenadas
        coo_r = run_cst_TCC2_3D(pop.v_ex_r,pop.v_in_r,dat,[dat.N1_r,dat.N2_r]); coo_r = coo_r[0:-1,:]*c_r
        coo_t = run_cst_TCC2_3D(pop.v_ex_t,pop.v_in_t,dat,[dat.N1_t,dat.N2_t]); coo_t = coo_t[0:-1,:]*c_t

        # Alguns dados a mais
        nb1 = dat.nb1
        if nb1 == 'L': # Tornar uniforme (ou próximo disso) a distribuição de seções ao longo da envergadura completa
            nb = dat.nb
            # Alterar a quantidade de seções em termos de uma concentração especificada
            if nb[1] == 0: # Usar o valor especificado 
                nb = nb[0]
            else: # Usar como uma concentração por metro e determinar o número de seções
                nb = int(np.floor((b*nb[0]-2)/2))
            
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
        
        sec_af_N = coo_r.shape[0] # Número de nós por seção (cada aerofólio)
        sec_N = 5 + nb1*2 + nb2*2 # Número de seções transversais
        far = b*2
            
        # Distribuição de torções geométricas
        if tw_m == 'L' and not isinstance(tw_t,str): # Definir a torção do meio em termos da torção da ponta
            tw_m = tw_t*b1/b
            tw_V = np.linspace(0,1,3+nb1+nb2)*tw_t
            tw_V1 = np.flip(tw_V[1:nb1+1])
            tw_V2 = np.flip(tw_V[2+nb1:-1])
        elif tw_t == 'L' and not isinstance(tw_m,str): # Definir a torção da ponta em termos da torção do meio
            tw_t = tw_m*b/b1
            tw_V = np.linspace(0,1,3+nb1+nb2)*tw_t
            tw_V1 = np.flip(tw_V[1:nb1+1])
            tw_V2 = np.flip(tw_V[2+nb1:-1])
        else: # Aplicar ambos os valores de torção especificados
            tw_V1 = np.linspace(0,1,2+nb1)*tw_m
            tw_V2 = np.linspace(tw_m/tw_t,1,2+nb2)*tw_t
            tw_V1 = np.flip(tw_V1[1:-1])
            tw_V2 = np.flip(tw_V2[1:-1])
        
        # (Em todos os casos acima, as torções intermediárias são definidas linearmente)

        # Rotacionar o perfil da ponta
        coo_t_R = airfoil_rotation(coo_t,tw_t)
        coo_t[:,0] += 1/4 # Consertar as coordenadas (coo_t é afetado pela função airfoil_rotation. não sei por quê)
        
        # Carregar as coordenadas do aerofólio do meio
        if 'L' in pop.v_ex_m or 'L' in pop.v_in_m:
            coo_m = airfoil_interpolation(coo_r,coo_t,b1/b,0,tw_t*b1/b)
        else:
            coo_m = run_cst_TCC2_3D(pop.v_ex_m,pop.v_in_m,dat,[dat.N1_m,dat.N2_m]); coo_m = coo_m[0:-1,:]*c_m

        # Obter nós ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        NODE = np.zeros((sec_af_N*(1+nb2),3))
        if sweep1 == 'Z': # Fazer com que o enflechamento da linha c/2 seja sempre zero
            dx_m = (c_r - c_m*cos(radians(tw_m)))/2 # Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
        else: # Usar enflechamento especificado pelo usuário
            dx_m = b1/2*tan(radians(sweep1))
        
        if sweep2 == 'Z': # Fazer com que o enflechamento da linha c/2 seja sempre zero
            dx_t = (c_m - c_t*cos(radians(tw_t)))/2 + dx_m # Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
            sweep2_angle = degrees(atan((c_t-c_m)/(b2/2)))
        else: # Usar enflechamento especificado pelo usuário
            dx_t = b2/2*tan(radians(sweep2)) + dx_m
            sweep2_angle = sweep2
        
    #    dx_t = (c_r - c_t*cosd(tw_t))/2; # Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
    #    dx_m = (c_r - c_m*cosd(tw_m))/2;
        NODE[0:sec_af_N,:] = np.hstack((np.array([coo_t_R[:,0]]).transpose()+dx_t,repmat(-b/2,sec_af_N,1),np.array([coo_t_R[:,1]]).transpose())) # Ponta esquerda da asa

        # Gerar seções intermediárias (ponta/meio)
        if nb2 > 0:
            # Criar struct que guarda seções de asa intermediárias
            # (isto será um auxílio devido à natureza simétrica a asa)
            wing_sec2 = structure()
            wing_sec2.coo = None; wing_sec2.dx = 0
            wing_sec2 = wing_sec2.repeat(nb2) # Inicializar
            op_vec2 = np.linspace(1,0,2+nb2); op_vec2 = op_vec2[1:-1] # Definição do formato da interpolação em função dos originais
            
            
    #        # usar este para sweep2 = 'Z'
    #        if sweep2 == 'Z'
    #dx_t = (c_m - c_t*cosd(tw_t))/2 + dx_m
    #            dx_t2 = (c_m - c_t*cosd(tw_t))/2; # Mesma coisa que dx_t, mas em relação à seção do meio
    ##            op_vec2b = op_vec2;
    #        else # usar este para sweep2 = número
    #            dx_t2 = b2/2*tand(sweep2) + dx_m;
    ##            op_vec2b = linspace(1,0,3+nb1+nb2);
    ##            op_vec2b = op_vec2b(2+nb1:end-1);
    #            
    #        end
        
            # Encontrar as coordenadas das seções (fazer interpolações)
            for i in range (len(op_vec2)):
                # Fazer interpolação e aplicar torção geométrica
                tw_sec = tw_V2[i]
                wing_sec2[i].coo = airfoil_interpolation(coo_m,coo_t,op_vec2[i],-op_vec2[i]*((b-b1)/2)-b1/2,tw_sec)
                # Aplicar o translado
                if sweep2 == 'Z':
                    wing_sec2[i].dx = op_vec2[i]*(dx_t-dx_m)+dx_m
                else:
                    wing_sec2[i].dx = op_vec2[i]*b2/2*tan(radians(sweep2_angle))+dx_m
                
                wing_sec2[i].coo[:,0] += wing_sec2[i].dx
                # Adicionar ao lado esquerdo da asa
                NODE[sec_af_N*(i+1):sec_af_N*(i+2),:] = wing_sec2[i].coo

        # Rotacionar perfil do meio
        coo_m_R = airfoil_rotation(coo_m,tw_m)
        coo_m[:,0] += 1/4

        # Adicionar coordenadas do perfil do meio
        NODE = np.vstack((NODE,np.hstack((np.array([coo_m_R[:,0]]).transpose()+dx_m,np.zeros((sec_af_N,1))-b1/2,np.array([coo_m_R[:,1]]).transpose()))))

        # Gerar seções intermediárias (meio/raiz)
        if nb1 > 0:
            temp = np.zeros((sec_af_N*(nb1),3))
            # Criar struct que guarda seções de asa intermediárias
            # (isto será um auxílio devido à natureza simétrica a asa)
            wing_sec1 = structure()
            wing_sec1.coo = None; wing_sec1.dx = 0
            wing_sec1 = wing_sec1.repeat(nb1) # Inicializar
            op_vec1 = np.linspace(1,0,2+nb1); op_vec1 = op_vec1[1:-1] # Definição do formato da interpolação em função dos originais
            
            # Encontrar as coordenadas das seções (fazer interpolações)
            for i in range(len(op_vec1)):
                # Fazer interpolação e aplicar torção geométrica
                tw_sec = tw_V1[i]
                wing_sec1[i].coo = airfoil_interpolation(coo_r,coo_m,op_vec1[i],-op_vec1[i]*(b1/2),tw_sec)
                # Aplicar o translado
                wing_sec1[i].dx = op_vec1[i]*dx_m
                wing_sec1[i].coo[:,0] += wing_sec1[i].dx
                # Adicionar ao lado esquerdo da asa
                temp[sec_af_N*i:sec_af_N*(i+1),:] = wing_sec1[i].coo
            
            NODE = np.vstack((NODE,temp))
        

        # Adicionar coordenadas do perfil da raiz
        NODE = np.vstack((NODE,np.hstack((np.array([coo_r[:,0]]).transpose(),np.zeros((sec_af_N,1)),np.array([coo_r[:,1]]).transpose()))))

        # Adicionar seções intermediárias ao lado direito da asa (raiz/meio)
        if nb1 > 0:
            temp = np.zeros((sec_af_N*nb1,3))
            k = 0
            for i in range(len(op_vec1)-1,-1,-1):
                wing_sec1[k].coo[:,1] = -wing_sec1[k].coo[:,1] # Inverter o sinal da coordenada y das seções intermediárias
                temp[sec_af_N*(i+1)-sec_af_N:sec_af_N*(i+1),:] = wing_sec1[k].coo
                k += 1
            
            NODE = np.vstack((NODE,temp))
        
        # Adicionar perfil do meio do lado direito da asa
        NODE = np.vstack((NODE,np.hstack((np.array([coo_m_R[:,0]]).transpose()+dx_m,np.zeros((sec_af_N,1))+b1/2,np.array([coo_m_R[:,1]]).transpose()))))

        # Adicionar seções intermediárias ao lado direito da asa (meio/ponta)
        if nb2 > 0:
        #    NODE = [NODE;zeros(sec_af_N*nb,3)];
            temp = np.zeros((sec_af_N*nb2,3))
            k = 0
            for i in range(len(op_vec2)-1,-1,-1):
                wing_sec2[k].coo[:,1] = -wing_sec2[k].coo[:,1] # Inverter o sinal da coordenada y das seções intermediárias
                temp[sec_af_N*(i+1)-sec_af_N:sec_af_N*(i+1),:] = wing_sec2[k].coo
                k += 1
            
            NODE = np.vstack((NODE,temp))

        # Adicionar ponta direita
        NODE = np.vstack((NODE,np.hstack((np.array([coo_t_R[:,0]]).transpose()+dx_t,np.zeros((sec_af_N,1))+b/2,np.array([coo_t_R[:,1]]).transpose()))))
        
        # Calcular área e corda aerodinâmica média
        S1 = (c_r + c_m)*b1/2; S2 = (c_m + c_t)*b2/2
        S = S1 + S2
        mac = (S1/b1 + S2/b2)/2 # obtido a partir das relações da razão de aspecto
        
        
    # Gerar painéis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Número de painéis: 
    panel_surf = np.zeros((sec_af_N*(sec_N-1),9)) # Inicializar painéis da superfície
    # Devido à natureza das geometrias, não haverão elementos triangulares 
    # Nota: painéis da superfície que se localizam no bordo de fuga não tem um dos
    # painéis adjacentes. Nesse caso, insere-se o valor 0

    # Montar os elementos da superfície
    v = np.arange(1,sec_N)
    k = 0
    for i in v:
        for j in range(sec_af_N-1):
            panel_surf[k,0:5] = np.hstack((1,k+1,sec_af_N+k+1,sec_af_N+k+2,k+2))
            k += 1
        
        panel_surf[k,0:5] = np.hstack((1,k+1,sec_af_N+k+1,sec_af_N*i+1,sec_af_N*(i-1)+1))
        k += 1

    # Adicionar numeração de painéis adjacentes
    panel_surf[0,5:] = np.hstack((sec_af_N+1,2,0,0)) # Primeiro painel
    panel_surf[-1,5:] = np.hstack((sec_af_N*(sec_N-2),sec_af_N*(sec_N-1)-1,0,0)) # Último painel
    for i in range(1,panel_surf.shape[0]-1):
        if i < sec_af_N-1: # Ponta esquerda da asa
            panel_surf[i,5:] = np.hstack((i,i+sec_af_N+1,i+2,0)); #print('caso1')
        elif i > panel_surf.shape[0]-sec_af_N: # Ponta direita da asa
            panel_surf[i,5:] = np.hstack((i,0,i+2,i-sec_af_N+1)); #print('caso2')
        elif i == sec_af_N-1: # Painel do intradorso da ponta esquerda (bordo de fuga)
            panel_surf[i,5:] = np.hstack((sec_af_N*2,0,0,sec_af_N-1)); #print('caso3')
        elif i == sec_af_N*(1+2*nb): # Painel do extradorso da ponta direita(bordo de fuga)
            panel_surf[i,5:] = np.hstack((0,0,sec_af_N*(1+nb*2)+2,sec_af_N*(1+nb)+1)); #print('caso4')
        else: # Todos os outros pontos
            panel_surf[i,5:] = np.hstack((i,i+sec_af_N+1,i+2,i-sec_af_N+1)); #print('caso5')

    if op == 1:
        plt.figure()
        plt.scatter(NODE[:,0],NODE[:,1])
        plt.axis('equal');plt.grid('true')

    # Adicionar os nós da trilha da asa
    # Cada um será posicionado diretamente atrás de sua respectiva seção de asa a uma distância far
    far_nodes = np.zeros((sec_N,3))
    for i in range(1,sec_N+1):
        far_nodes[i-1,:] = np.hstack((NODE[sec_af_N*i-sec_af_N,0]+far,NODE[sec_af_N*i-sec_af_N,1:3]))

    A = NODE.shape[0]
    NODE = np.vstack((NODE,far_nodes))

    # Montar os painéis da trilha
    panel_far = np.zeros((sec_N-1,9))
    k = 1
    for i in range(sec_N-1):
        panel_far[i,:] = np.hstack((10,k,A+i+1,A+i+2,k+sec_af_N,k,k+sec_af_N-1,0,0))
        k += sec_af_N

    PANEL = np.vstack((panel_surf,panel_far))
    
    # Inserir as novas informações ao struct original
    pop.NODE = NODE
    pop.PANEL = PANEL
    pop.sec_N = sec_N
    pop.S = S
    pop.mac = mac
    
    return pop