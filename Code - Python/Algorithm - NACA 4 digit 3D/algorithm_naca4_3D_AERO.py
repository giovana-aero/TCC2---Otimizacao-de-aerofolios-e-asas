## Algoritmo genético
# Otimização de asas trapezoidais simples ou duplas, com perfis NACA 4 dígitos

# Nota pra posterioridade: grandes quantidades de pontos nos aerofólios (np)
# tendem a fazer o apame travar de vez em quando




# A FAZER: revisar a configuração da malha (u1, w1, essas coisas)
#       é preciso traduzir a fução de interpolação de funções devido à opção 'L' do aerofólio do meio
          


# CONSERTAR ISTO

#  - Dados dos aerofólios - 
# Raiz: NACA 3.04.010.0
# Meio: NACA 3.03.016.0
# Ponta: NACA 2.03.012.0
# Torção geométrica na ponta tw_t = -4.800000000000001
#  - Dados aerodinâmicos - 
# Condição de voo 0: 0.047720, 0.001160, 41.006140, -0.054310
# Condição de voo 1: 0.000000, 0.000000, 0.000000, 0.000000
# Pontuação: 1.454400






# Pacotes necessários
import numpy as np
from ypstruct import structure
import random
# import os
import matplotlib.pyplot as plt
from time import time
import openvsp as vsp
# Funções utilizadas
from error_check_naca4_3D import error_check_naca4_3D
from run_apame_mesher_naca4 import run_apame_mesher_naca4
from run_apame import run_apame
# from run_openvsp_naca4 import run_openvsp_naca4
from run_openvsp_naca4_VLM import run_openvsp_naca4_VLM
from run_openvsp_naca4_VLM_v2 import run_openvsp_naca4_VLM_v2
# from run_vspaero import run_vspaero
from run_vspaero_VLM import run_vspaero_VLM
from fitness_naca4_3D import fitness_naca4_3D
# from make_vector import make_vector
from selection_crossover import selection_crossover
from plot_planform import plot_planform


start = time()


# Parâmetros do algoritmo
dat = structure()
dat.N = 100                            # Número de indivíduos na população
dat.mu = 0.05                           # Probabilidade de mutação (definida entre zero e um)
dat.iter = 5                           # Número de iterações
dat.elite = 1                         # Aplicar elitismo?
dat.subs = 1                          # Substituir asas sem resultados? (ver ainda se isto será necessário)

# Parâmetros da geometria: planta da asa
dat.planf_op = 0.5 # Proporção de asas trapezoidais simples e bitrapezoidais (0->todas trapezoidais simples, 1->todas bitrapezoidais)
dat.b_ext_in = [10,15] # Envergadura completa [m] (extensão inicial)
dat.b_step = 0.5 # Valor do passo para definição das extensões de envergadura [m]
dat.b1_ext_in = [8,14,0.5] # Envergadura da raiz ao meio [m] (asas bitrapezoidais apenas) (valor mínimo,valor máximo,separação mínima da ponta da asa (considerando apenas uma metade)) (extensão inicial)
dat.b1_step = 0.5 # Valor do passo para definição das extensões da envergadura da primeira seção [m]
dat.c_r_ext_in = [1,2] # Corda da raiz [m] (extensão inicial)
dat.c_r_step = 0.1 # Valor do passo para definição das extensões da corda da raiz [m]
dat.c_m_ext_in = [0.5,2] # Corda do meio [m] (asas bitrapezoidais apenas) (extensão inicial) (a opção 'L' força o formato trapezoidal simples)
dat.c_m_step = 0.1 # Valor do passo para definição das extensões da corda do meio [m]
dat.c_t_ext_in = [0.5,2] # Corda da ponta [m] (extensão inicial)
dat.c_t_step = 0.1 # Valor do passo para definição das extensões da corda da raiz [m]
dat.sweep_ext_in = [0,15] # Enflechamento de asas trapezoidais simples (extensão inicial) (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep_step = 1 # Valor do passo para definição das extensões do enflechamento
dat.sweep1_ext_in = [0,15] # Enflechamento da primeira seção de asas trapezoidais duplas(extensão inicial) (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep1_step = 1 # Valor do passo para definição das extensões do enflechamento
dat.sweep2_ext_in = [0,15] # Enflechamento da segunda seção de asas trapezoidais duplas(extensão inicial) (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep2_step = 1 # Valor do passo para definição das extensões do enflechamento
# Parâmetros da geometria: aerofólio da raiz
dat.m_ext_r = [0,4] # Curvatura máxima
dat.p_ext_r = [1,4] # Local da curvatura máxima
dat.t_ext_r = [10,20] # Espessura máxima 
# Parâmetros da geometria: aerofólio do meio (asas bitrapezoidais apenas) 
dat.m_ext_m = [0,4] # (a opção 'L' aqui gera o perfil do meio linearmente em função dos perfis da raiz e da ponta)
dat.p_ext_m = [1,4]
dat.t_ext_m = [10,20]
# Parâmetros da geometria: aerofólio da ponta
dat.m_ext_t = [0,4]
dat.p_ext_t = [1,4]
dat.t_ext_t = [10,20]
dat.tw_t_ext_in = [-5,0] # Torção geométrica [°]
dat.tw_t_step = 0.25 # Valor do passo para definição das extensões da torção geométrica na ponta [°]

# Parâmetros da malha
dat.np = 30 # Número de pontos na geração de ordenadas nos aerofólios
dat.np_op = 1 # 1 -> cosspace, 0 -> cosspace_half
dat.nb = [2,1] # Número de seções intermediárias (raiz/ponta) [número de seções,0] ou [concentração por metro,1]
dat.nb1 = 'L' # Número de seções intermediárias (raiz/meio) (asas bitrapezoidais apenas) (opção 'L' faz com que nb1 e nb2 sejam uniformemente determinados ao longo da envergadura)
dat.nb2 = 0 # Número de seções intermediárias (meio/ponta) (asas bitrapezoidais apenas)
# dat.num_U1 = 10
# dat.num_w1 = 10

# Parâmetros das simulações
dat.method = 2                        # Método da simulação (1->painéis 3D no APAME, 2->VLM no VSPAERO)
dat.cases = 1                          # Número de condições de voo a serem analisadas
dat.v_ref = [100,100,100] # Velocidades de referência [m/s] 
dat.rho = [1.225,1.225,1.225] # Densidades do ar [kg/m^3] 
dat.p_atm = [101325,101325,101325] # Pressões do ar [Pa] (irrelevante neste algoritmo)
dat.mach = [0,0.,0.2] # Números de Mach
dat.reynolds = [1e6,1e6,1e6]           # Valores dos números de Reynolds para as simulações (irrelevante neste algoritmo))
dat.aoa = [5,2,4]                     # Ângulos de ataque
dat.karman_tsien = ['N','N','N']    # Correção de compressibilidade de Karman-Tsien
dat.wake_iters = [1,1,1]            # Números de iterações no VSPAERO
dat.coeff_op = np.array((['!','!','^','!'],       # Uma linha para cada condição de voo
                         ['!','!','!','!'],
                         ['!','!','!','!']))
dat.coeff_val = np.array(([0.05,7e-3,90,-1e-1],
                          [0.5,0,0,-0.08],
                          [0,0,0,-0.08]))
dat.coeff_F = np.array(([1,1,1,1],
                        [1,1,1,1],
                        [1,1,1,1]))
# [CL CD L/D CM] Definição de cada linha da matriz dat.coeff_op
# '!' -> não usar como função objetiva
# '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD)
# 'c' -> buscar valor constante de coeficiente de momento (arbitrário)
# 'k' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val(1,4))
# 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
# em sua respectiva casa na matriz dat.coeff_val
# 'q' -> procurar por um valor específico de força de sustentação, força de arrasto
# ou momento de arfagem (CL, CD e CM). Nesse caso, definir o valor em sua casa na matriz dat.coeff_val
# '#' -> procurar por um valor máximo (L) ou mínimo (D)
# A matriz dat.coeff_F dá os pesos de cada função objetiva

# Checagem de erros
dat = error_check_naca4_3D(dat)

# Template dos structs
# Forma da planta
empty = structure()
empty.type = None
empty.b = None
empty.b1 = None 
empty.c_r = None 
empty.c_m = None 
empty.c_t = None
empty.mac = None
empty.S = None
empty.sweep = None
empty.sweep1 = None
empty.sweep2 = None
# Aerofólios (todos serão definidos com vetores [m,p,t])
empty.af_r = np.zeros(3)
empty.af_m = np.zeros(3)
empty.af_t = np.zeros(3)
empty.tw_m = None
empty.tw_t = None
# Dados da malha
#empty.NODE = None
#empty.PANEL = None
#empty.sec_N = None
# Dados aerodinâmicos e pontuação
empty.aero = None # Terá o mesmo formato que a matriz coeff_op
empty.score = 0

# Inicializar
pop = empty.repeat(dat.N)
chi = empty.repeat(dat.N)

# Redefinir as extensões dos valores das variáveis a partir das extensões iniciais
# (precisão dos valores podem ser alteradas aqui)
dat.b_ext = np.arange(dat.b_ext_in[0],dat.b_ext_in[1]+dat.b_step,dat.b_step)
dat.c_r_ext = np.arange(dat.c_r_ext_in[0],dat.c_r_ext_in[1]+dat.c_r_step,dat.c_r_step)
if not 'Z' in dat.sweep_ext_in: dat.sweep_ext = np.arange(dat.sweep_ext_in[0],dat.sweep_ext_in[1]+dat.sweep_step,dat.sweep_step)
if not 'Z' in dat.sweep1_ext_in: dat.sweep1_ext = np.arange(dat.sweep1_ext_in[0],dat.sweep1_ext_in[1]+dat.sweep1_step,dat.sweep1_step)
if not 'Z' in dat.sweep2_ext_in: dat.sweep2_ext = np.arange(dat.sweep2_ext_in[0],dat.sweep2_ext_in[1]+dat.sweep2_step,dat.sweep2_step)
dat.tw_t_ext = np.arange(dat.tw_t_ext_in[0],dat.tw_t_ext_in[1],dat.tw_t_step)
# As demais extensões de valores serão definidos para cada indivíduo a seguir 
# de modo a cumprir os seguintes requisitos:
# b1 < b
# c_r >= c_m >= c_t
# ATENÇÃO: PÔR ISTO NA CHECAGEM DE ERROS DEPOIS

# Se usar o VSPAERO, gerar a asa preliminar (terá planta trapezoidal simples)
if dat.method == 2:
    # Iniciar modelo
    vsp.ClearVSPModel()
    vsp.VSPRenew()

    # Adicionar asa
    wid = vsp.AddGeom( "WING", "" )
    
    # Adicionar mais uma seção
    # vsp.InsertXSec( wid, 1, vsp.XS_FOUR_SERIES )

    # Mudar drivers
    # isto funciona da mesma forma que no gui: o projetista determina três parâmetros
    # entre os 8 possíveis pra manipular a planta da asa
    vsp.SetDriverGroup( wid, 1, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    # vsp.SetDriverGroup( wid, 2, vsp.SPAN_WSECT_DRIVER, vsp.ROOTC_WSECT_DRIVER, vsp.TIPC_WSECT_DRIVER )
    
    vsp.Update()
    # last_planf = 0
    last_planf = 0

# Gerar população inicial
print('<< Geração da população inicial >>')
for i in range(dat.N):
    print('Indivíduo ' + str(i+1))
    
    # Estabelecer tipo da planta
    # 0 -> trapezoidal simples, 1 -> bitrapezoidal
    pop[i].type = int(random.random() <= dat.planf_op)
    
    # Gerar forma da planta
    pop[i].b = dat.b_ext[random.randint(0,len(dat.b_ext)-1)]
    pop[i].c_r = dat.c_r_ext[random.randint(0,len(dat.c_r_ext)-1)]
    if dat.c_t_ext_in[1] > pop[i].c_r:
        # Caso, nesta asa, o valor máximo da extensão da corda da ponta seja 
        # maior do que a corda da raiz, usar a corda da raiz atual como 
        # o valor máximo da extensão
        dat.c_t_ext = np.arange(dat.c_t_ext_in[0],pop[i].c_r,dat.c_t_step)  #np.arange(dat.c_t_ext_in[0],pop[i].c_r+dat.c_t_step,dat.c_t_step)
    else:
        # Caso contrário, usar os valores definidos pelo usuário em dat.c_t_ext_in
        dat.c_t_ext = np.arange(dat.c_t_ext_in[0],dat.c_t_ext_in[1]+dat.c_t_step,dat.c_t_step)
    pop[i].c_t = dat.c_t_ext[random.randint(0,len(dat.c_t_ext)-1)]
    pop[i].tw_t = dat.tw_t_ext[random.randint(0,len(dat.tw_t_ext)-1)]
    
    # Gerar aerofólio da raiz
    pop[i].af_r[0] = random.randint(dat.m_ext_r[0],dat.m_ext_r[1]) # Curvatura máxima
    if pop[i].af_r[0] == 0: # Perfis simétricos
        pop[i].af_r[1] = 0
    else: # Perfis assimétricos
        pop[i].af_r[1] = random.randint(dat.p_ext_r[0],dat.p_ext_r[1])
    pop[i].af_r[2] = random.randint(dat.t_ext_r[0],dat.t_ext_r[1]) # Espessura máxima
    
    # Gerar aerofólio da ponta
    pop[i].af_t[0] = random.randint(dat.m_ext_t[0],dat.m_ext_t[1]) # Curvatura máxima
    if pop[i].af_t[0] == 0: # Perfis simétricos
        pop[i].af_t[1] = 0
    else: # Perfis assimétricos
        pop[i].af_t[1] = random.randint(dat.p_ext_t[0],dat.p_ext_t[1])
    pop[i].af_t[2] = random.randint(dat.t_ext_t[0],dat.t_ext_t[1]) # Espessura máxima
    
    # Dados adicionais para asas bitrapezoidais
    if pop[i].type == 1:
        # Mais dados da planta
        if 'L' in dat.c_m_ext_in:
            # Caso c_m seja definido como 'L', seu valor real será atribúido pela
            # função run_apame_mesher_naca4
            pop[i].c_m = 'L'
        else:
            # Caso contrário, é realizado o processo abaixo
            dat.c_m_ext = [0,0]
            if dat.c_m_ext_in[0] < pop[i].c_t:
                # Caso o valor mínimo da extensão da corda do meio seja menor do que
                # a corda da ponta da asa atual, usar o valor da corda da ponta da
                # asa como o valor mínimo da extensão
                dat.c_m_ext[0] = pop[i].c_t
            else:
                # Caso contrário, usar o valor original
                dat.c_m_ext[0] = dat.c_m_ext_in[0]
            
            if dat.c_m_ext_in[1] > pop[i].c_r:
                # Caso o valor mpaximo da extensão da corda do meio seja maior do que
                # a corda da raiz da asa atual, usar o valor da corda da raiz da
                # asa como o valor mínimo da extensão
                dat.c_m_ext[1] = pop[i].c_r
            else:
                # Caso contrário, usar o valor original
                dat.c_m_ext[1] = dat.c_m_ext_in[1]
            
            pop[i].c_m = dat.c_m_ext[random.randint(0,len(dat.c_m_ext)-1)]
        
        
        if dat.b1_ext_in[1] > pop[i].b-dat.b1_ext_in[2]*2:
            # Se, nesta asa, o valor máximo da extensão de b1 for maior do que
            # o valor máximo permitido pelo requisito de separação, utilizar
            # a envergadura da asa atual menos a separação mínima como valor
            # máximo da extensão
            # O valor da separação mínima é utilizado para garantir que b1 ~= b
            dat.b1_ext = np.arange(dat.b1_ext_in[0],pop[i].b-dat.b1_ext_in[2]*2+dat.b1_step,dat.b1_step)
        else:
            # Caso contrário, utilizar a extensão original definida em dat.b1_ext_in
            dat.b1_ext = np.arange(dat.b1_ext_in[0],dat.b1_ext_in[1]+dat.b1_step,dat.b1_step)
        
        pop[i].b1 = dat.b1_ext[random.randint(0,len(dat.b1_ext)-1)]

        # Aerofólio do meio
        pop[i].tw_m = 'L' # Essa sempre será a configuração deste algoritmo, mas isso pode ser alterado (com as devidas alterações no resto do código)
        pop[i].af_m[0] = random.randint(dat.m_ext_m[0],dat.m_ext_m[1]) # Curvatura máxima
        if pop[i].af_m[0] == 0: # Perfis simétricos
            pop[i].af_m[1] = 0
        else: # Perfis assimétricos
            pop[i].af_m[1] = random.randint(dat.p_ext_m[0],dat.p_ext_m[1])
        pop[i].af_m[2] = random.randint(dat.t_ext_m[0],dat.t_ext_m[1]) # Espessura máxima
        
        
        # Enflechamento da primeira seção
        if not 'Z' in dat.sweep1_ext_in:
            pop[i].sweep1 = dat.sweep1_ext[random.randint(0,len(dat.sweep1_ext)-1)]
        else:
            pop[i].sweep1 = 'Z'
        
        # Enflechamento da segunda seção
        if not 'Z' in dat.sweep2_ext_in:
            pop[i].sweep2 = dat.sweep2_ext[random.randint(0,len(dat.sweep2_ext)-1)]
        else:
            pop[i].sweep2 = 'Z'
        
    else: # Estabelecer enflechamento da asa trapezoidal simples
        if not 'Z' in dat.sweep_ext_in:
            pop[i].sweep = dat.sweep_ext[random.randint(0,len(dat.sweep_ext)-1)]
        else:
            pop[i].sweep = 'Z'
            
# Gerar struct que guarda a melhor asa de cada geração
archive = empty.repeat(dat.iter)

## Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for loop in range(dat.iter):
    print('Iteração ' + str(loop+1))
    
    # Simular as asas e obter dados
    select = np.ones((1,dat.N))
    print('<< Simulação das asas >>')
    if dat.method == 1: # Painéis 3D
        for i in range(dat.N):
            pop[i] = run_apame_mesher_naca4(pop[i],dat,1) # Obter malha do apame
            pop[i].aero = run_apame(pop[i],dat) # Fazer simulação
            
            # Marcar indivíduos que não tenham convergido na simulação
            if isinstance(pop[i].aero,str) and pop[i].aero == 'n':
                select[[0],[i]] = 0
            
    elif dat.method == 2: # VLM
        for i in range(dat.N):
            print('Indivíduo ' + str(i+1))
            # pop[i] = run_openvsp_naca4_VLM(pop[i],dat) # Obter geometria OpenVSP            
            pop[i],last_planf,flag = run_openvsp_naca4_VLM_v2(pop[i],dat,wid,last_planf,i) # Obter geometria OpenVSP            
            # print('simulação')
            if flag == 1: 
                pop[i].aero = run_vspaero_VLM(pop[i],dat) # Simulação
            else:
                pop[i].aero = 'n'
            
            # Marcar indivíduos que não tenham convergido na simulação
            if isinstance(pop[i].aero,str) and pop[i].aero == 'n':
                select[[0],[i]] = 0
    
    # Encontrar indivíduos problemáticos para já ignorá-los durante a atribuição de pontuação
    # select = np.expand_dims(np.argwhere(select==0)[:,1],axis=0)
    select = np.argwhere(select==0)[:,1]
    select2 = np.arange(0,dat.N)
    select2 = np.delete(select2,select)

    if len(select) == dat.N:
        raise TypeError('Nenhum aerofólio convergiu nas simulações')
    
    # Atribuir pontuações (fitnesses)
    pop = fitness_naca4_3D(pop,dat,select2)

    # Pôr todas as pontuações em um vetor
    weights = np.zeros(dat.N)
    for i in range(dat.N):
        weights[i] = pop[i].score
    
    # Guardar o melhor perfil de cada iteração
    pos = np.argmax(weights)
    archive[loop] = pop[pos].deepcopy()
    
    # Mostrar a melhor asa
    plt.figure()
    plot_planform(pop[pos],'k')
    plt.axis('equal');plt.grid('true')
    plt.title('Iteração ' + str(loop+1))
    
    # Parar o código aqui na última iteração, já que nesse cenário o resto
    # do código é inútil
    if loop == dat.iter-1: break

    # Substituir indivíduos com pontuação nula ou negativa por aqueles com as
    # pontuações mais altas
    if dat.subs == 1:
        # Agora o vetor select também inclui indivíduos que convergiram na
        # simulação, mas que são extremamente inaptos (pontuação negativa)
        select = np.argwhere(weights <= 0)
        if select.size != 0: # Se select não estiver vazio

            if len(select) <= len(select2): # Se o número de aerofólios com aero = 'n' for menor ou igual que o número de aerofólios com pontuação
                # Preencher o vetor ind com os índices dos indivíduos de maior pontuação
                ind = np.zeros(len(select))
                temp = weights
                for i in range(len(ind)):
                    ind[i] = np.argmax(temp)
                    temp[int(ind[i])] = -np.inf
                # Substituir os indivíduos
                k = 0
                for i in select[0,:]:
                    pop[i] = pop[int(ind[k])].deepcopy()
                    k += 1
            else: 
                # Preencher o vetor ind com os índices dos indivíduos de maior pontuação (que acabam sendo todos aqueles apontados por select2)
                ind = select2
                # Substituir os indivíduos
                k = 0
                for i in select[:,0]:
                    pop[i] = pop[ind[k]].deepcopy()
                    k += 1
                    if k == len(ind): k = 0

    # Escolher membros da população pra reprodução
    c = 0
    print('<< Reprodução >>')
    for f in range (int(dat.N/2)):
        print('{:.2f}% completo'.format((f+1)/(dat.N/2)*100))
		
        # Isto seleciona dois pais por meio de uma seleção via roleta
        # (indivíduos com pesos maiores têm  mais chance de serem selecionados)
        par = [0,0] # Vetor que indica a numeração dos pais escolhidos
        par[0] = selection_crossover(weights)
        par[1] = selection_crossover(weights)
        
        # Crossover, gerando dois filhos de cada dois pais
        # Genes a serem trocados
        # Envergadura b
        # Envergadura b1 (asas bitrapezoidais apenas)
        # Corda da raiz
        # Corda do meio (asas bitrapezoidais apenas)
        # Corda da ponta
        # - aerofólio da raiz
        # - aerofólio do meio (asas bitrapezoidais apenas)
        # - aerofólio da ponta
        # Torção geométrica na ponta
        
        if pop[par[0]].type != pop[par[1]].type: # Caso as asas sejam de tipos diferentes
            op = random.randint(1,2)
        else: # Caso as asas sejam de tipos iguais
            op = 1        

        if op == 1: # Fazer crossover normalmente (trocar genes gerais)
            # Escrever um vetor de opções:
            # [b,b1,c_r,c_m,c_t,airfoil_r,airfoil_m,airfoil_t,tw_t]
            cross_op_v = np.random.randint(2,size=12)
            
            chi[c] = pop[par[0]]
            chi[c+1] = pop[par[1]]
            
            # Envergadura b
            if cross_op_v[0] == 1:
                chi[c].b = pop[par[1]].b               
                chi[c+1].b = pop[par[0]].b
        
            # Envergadura da primeira seção b1 [asas bitrapezoidais apenas] [deve cumprir o requisito de separação mínima]
            if cross_op_v[1] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1: #&& pop[par[1]].b1 <= pop[par[2]].b - dat.b1_ext_in[3]*2 && pop[par[2]].b1 <= pop[par[1]].b - dat.b1_ext_in[3]*2
                chi[c].b1 = pop[par[1]].b1
                chi[c+1].b1 = pop[par[0]].b1
            
            # Corda da raiz c_r [deve cumprir o requisito c_r >= c_t]
#            if cross_op_v[2] == 1 and pop[par[0]].c_r >= pop[par[1]].c_t and pop[par[1]].c_r >= pop[par[0]].c_t:
            if cross_op_v[2] == 1 and chi[c].c_r >= chi[c+1].c_t and chi[c+1].c_r >= chi[c].c_t:
                chi[c].c_r = pop[par[1]].c_r
                chi[c+1].c_r = pop[par[0]].c_r
            
            # Corda da ponta c_t [deve cumprir o requisito c_r >= c_t]
#            if cross_op_v[3] == 1 and pop[par[0]].c_r >= pop[par[1]].c_t and pop[par[1]].c_r >= pop[par[0]].c_t:
            if cross_op_v[3] == 1 and chi[c].c_r >= chi[c+1].c_t and chi[c+1].c_r >= chi[c].c_t:
                chi[c].c_t = pop[par[1]].c_t
                chi[c+1].c_t = pop[par[0]].c_t
                
            # Corda do meio c_m [asas bitrapezoidais apenas] [deve cumprir o requisito c_r >= c_m >= c_t] [ignorar caso a opção 'L' seja aplicada a c_m]
            if cross_op_v[4] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1 and not 'L' in dat.c_m_ext_in:
                chi[c].c_m = pop[par[1]].c_m
                chi[c+1].c_m = pop[par[0]].c_m
            
            # Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi[c].type == 1:
                
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[c].b1 > chi[c].b-dat.b1_ext_in[2]*2:
                    chi[c].b1 = chi[c].b-dat.b1_ext_in[2]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[c].c_m > chi[c].c_r and not 'L' in dat.c_m_ext_in:
                    chi[c].c_m = chi[c].c_r
            
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[c].c_m < chi[c].c_t and not 'L' in dat.c_m_ext_in:
                    chi[c].c_m = chi[c].c_t
                
            if chi[c+1].type == 1:
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[c+1].b1 > chi[c+1].b-dat.b1_ext_in[2]*2:
                    chi[c+1].b1 = chi[c+1].b-dat.b1_ext_in[2]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[c+1].c_m > chi[c+1].c_r and not 'L' in dat.c_m_ext_in:
                    chi[c+1].c_m = chi[c+1].c_r
                
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[c+1].c_m < chi[c+1].c_t and not 'L' in dat.c_m_ext_in:
                    chi[c+1].c_m = chi[c+1].c_t
            
            # Aerofólio da raiz
            if cross_op_v[5] == 1:
            
                af1 = chi[c].af_r
                af2 = chi[c+1].af_r
                
                if af1[0] == 0 or af2[0] == 0: # Caso um deles seja simétrico
                    n = random.randint(1,2)    
                else:
                    n = random.randint(1,4)
                
                if n == 1: # Trocar os perfis completamente
                    chi[c].af_r = af2
                    chi[c+1].af_r = af1
                
                if n == 2: # Trocar a espessura máxima
                    temp1 = np.hstack(([af1[0:2],af2[-1]]))
                    temp2 = np.hstack(([af2[0:2],af1[-1]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_r = temp1
                        chi[c+1].af_r = temp2
                    else:
                        chi[c].af_r = temp2
                        chi[c+1].af_r = temp1
                
                if n == 3: # Trocar o local da curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1],af1[2]]))
                    temp2 = np.hstack(([af2[0],af1[1],af2[2]]))
                    
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_r = temp1
                        chi[c+1].af_r = temp2
                    else:
                        chi[c].af_r = temp2
                        chi[c+1].af_r = temp1
                
                if n == 4: # Trocar a curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1:]]))
                    temp2 = np.hstack(([af2[0],af1[1:]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_r = temp1
                        chi[c+1].af_r = temp2
                    else:
                        chi[c].af_r = temp2
                        chi[c+1].af_r = temp1
            
            # Aerofólio do meio (asas bitrapezoidais apenas) (ignorar caso a opção 'L' seja aplicada ao aerofólio do meio)
            if cross_op_v[6] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1 and 'L' != dat.m_ext_m and 'L' != dat.p_ext_m and 'L' != dat.t_ext_m:
                af1 = chi[c].af_m
                af2 = chi[c+1].af_m
                
                if af1[0] == 0 or af2[0] == 0: # Caso um deles seja simétrico
                    n = random.randint(1,2)    
                else:
                    n = random.randint(1,4)
                
                if n == 1: # Trocar os perfis completamente
                    chi[c].af_m = af2
                    chi[c+1].af_m = af1
                
                if n == 2: # Trocar a espessura máxima
                    temp1 = np.hstack(([af1[0:2],af2[-1]]))
                    temp2 = np.hstack(([af2[0:2],af1[-1]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_m = temp1
                        chi[c+1].af_m = temp2
                    else:
                        chi[c].af_m = temp2
                        chi[c+1].af_m = temp1
                
                if n == 3: # Trocar o local da curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1],af1[2]]))
                    temp2 = np.hstack(([af2[0],af1[1],af2[2]]))
                    
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_m = temp1
                        chi[c+1].af_m = temp2
                    else:
                        chi[c].af_m = temp2
                        chi[c+1].af_m = temp1
                
                if n == 4: # Trocar a curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1:]]))
                    temp2 = np.hstack(([af2[0],af1[1:]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_m = temp1
                        chi[c+1].af_m = temp2
                    else:
                        chi[c].af_m = temp2
                        chi[c+1].af_m = temp1
                        
            # Aerofólio da ponta
            if cross_op_v[7] == 1 :
                af1 = chi[c].af_t
                af2 = chi[c+1].af_t
                
                if af1[0] == 0 or af2[0] == 0: # Caso um deles seja simétrico
                    n = random.randint(1,2)    
                else:
                    n = random.randint(1,4)
                
                if n == 1: # Trocar os perfis completamente
                    chi[c].af_t = af2
                    chi[c+1].af_t = af1
                
                if n == 2: # Trocar a espessura máxima
                    temp1 = np.hstack(([af1[0:2],af2[-1]]))
                    temp2 = np.hstack(([af2[0:2],af1[-1]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_t = temp1
                        chi[c+1].af_t = temp2
                    else:
                        chi[c].af_t = temp2
                        chi[c+1].af_t = temp1
                
                if n == 3: # Trocar o local da curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1],af1[2]]))
                    temp2 = np.hstack(([af2[0],af1[1],af2[2]]))
                    
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_t = temp1
                        chi[c+1].af_t = temp2
                    else:
                        chi[c].af_t = temp2
                        chi[c+1].af_t = temp1
                
                if n == 4: # Trocar a curvatura máxima
                    temp1 = np.hstack(([af1[0],af2[1:]]))
                    temp2 = np.hstack(([af2[0],af1[1:]]))
                    if random.randint(0,1) == 1: # Decidir pra qual asa vai cada um dos aerofólios novos
                        chi[c].af_t = temp1
                        chi[c+1].af_t = temp2
                    else:
                        chi[c].af_t = temp2
                        chi[c+1].af_t = temp1
            
            # Torção geométrica na ponta
            if cross_op_v[8] == 1: 
                chi[c].tw_t = pop[par[1]].tw_t
                chi[c+1].tw_t = pop[par[0]].tw_t
            
            # Enflechamento total (asas trapezoidais simples)
            if cross_op_v[9] == 1 and pop[par[0]].type == 0 and pop[par[1]].type == 0:
                chi[c].sweep = pop[par[1]].sweep
                chi[c+1].sweep = pop[par[0]].sweep
        
            # Enflechamento da primeira seção (asas bitrapezoidais)
            if cross_op_v[10] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1:
                chi[c].sweep1 = pop[par[1]].sweep1
                chi[c+1].sweep1 = pop[par[0]].sweep1
            
            # Enflechamento da segunda seção (asas bitrapezoidais)
            if cross_op_v[11] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1:
                chi[c].sweep2 = pop[par[1]].sweep2
                chi[c+1].sweep2 = pop[par[0]].sweep2
        
        else: # Transformar uma asa trapezoidal simples em trapezoidal dupla e vice-versa
            # (adiciona-se/retira-se o aerofólio do meio)
            if pop[par[0]].type == 0: # Se a primeira do par for trapezoidal simples
                chi[c] = pop[par[0]]
                chi[c].type = 1
                chi[c].b1 = pop[par[1]].b1
                chi[c].c_m = pop[par[1]].c_m
                chi[c].af_m = pop[par[1]].af_m
                chi[c].tw_m = pop[par[1]].tw_m
                chi[c].sweep1 = pop[par[1]].sweep1
                chi[c].sweep2 = pop[par[1]].sweep2
                
                chi[c+1] = pop[par[1]]
                chi[c+1].type = 0
                chi[c+1].b1 = None
                chi[c+1].c_m = None
                chi[c+1].af_m = np.zeros((1,3))
                chi[c+1].tw_m = None
                chi[c+1].sweep = pop[par[0]].sweep
                
            else: # Se a primeira do par for trapezoidal dupla
                chi[c] = pop[par[1]]
                chi[c].type = 1
                chi[c].b1 = pop[par[0]].b1
                chi[c].c_m = pop[par[0]].c_m
                chi[c].af_m = pop[par[0]].af_m
                chi[c].tw_m = pop[par[0]].tw_m
                chi[c].sweep1 = pop[par[0]].sweep1
                chi[c].sweep2 = pop[par[0]].sweep2
                
                chi[c+1] = pop[par[0]]
                chi[c+1].type = 0
                chi[c+1].b1 = None
                chi[c+1].c_m = None
                chi[c+1].af_m = np.zeros(3)
                chi[c+1].tw_m = None
                chi[c+1].sweep = pop[par[1]].sweep
             
            # Checar o requisito das envergaduras
            if chi[c].b1 > chi[c].b - dat.b1_ext_in[2]*2:
                # Se não cumprir o requisito, atribuir um novo valor
                if chi[c].b - dat.b1_ext_in[2]*2 >= dat.b1_ext_in[1]:
                    dat.b1_ext = np.arange(dat.b1_ext_in[0],dat.b1_ext_in[1]+dat.b1_step,dat.b1_step)
                else:
                    dat.b1_ext = np.arange(dat.b1_ext_in[0],chi[c].b-dat.b1_ext_in[2]*2+dat.b1_step,dat.b1_step)
                    
                chi[c].b1 = dat.b1_ext[random.randint(0,len(dat.b1_ext)-1)]
            
            # Checar o requisito dos comprimentos de corda
            if chi[c].c_t > chi[c].c_m and not 'L' in dat.c_m_ext_in or chi[c].c_m > chi[c].c_r and not 'L' in dat.c_m_ext_in:
                # Se não cumprir o requisito, atribuir um novo valor
                dat.c_m_ext = [0,0]
                if chi[c].c_t > dat.c_m_ext_in[0]:
                    dat.c_m_ext[0] = chi[c].c_t
                else:
                    dat.c_m_ext[0] = dat.c_m_ext_in[0]
                
                if chi[c].c_r < dat.c_m_ext_in[1]:
                    dat.c_m_ext[1] = chi[c].c_r
                else:
                    dat.c_m_ext[1] = dat.c_m_ext_in[1]
                    
                
                if np.arange(dat.c_m_ext[0],dat.c_m_ext[1],dat.c_m_step).size == 0:
                    chi[c].c_m = dat.c_m_ext[0]
                else:
                    dat.c_m_ext = np.arange(dat.c_m_ext[0],dat.c_m_ext[1],dat.c_m_step)  # np.arange(dat.c_m_ext[0],dat.c_m_ext[1]+dat.c_m_step,dat.c_m_step)
                    chi[c].c_m = dat.c_m_ext[random.randint(0,len(dat.c_m_ext)-1)]
                
            if chi[c+1].c_r < chi[c+1].c_t:
                chi[c+1].c_r = chi[c+1].c_t
                
                # dat.c_m_ext = np.arange(dat.c_m_ext[0],dat.c_m_ext[1],dat.c_m_step)   #np.arange(dat.c_m_ext[0],dat.c_m_ext[1]+dat.c_m_step,dat.c_m_step)
                # chi[c].c_m = dat.c_m_ext[random.randint(0,len(dat.c_m_ext)-1)]
        
        # Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
        if chi[c].type == 1:
            if chi[c].b1 > chi[c].b-2*dat.b1_ext_in[2] or chi[c].c_t > chi[c].c_m or chi[c].c_m > chi[c].c_r:
                raise TypeError('Problema na geometria da planta (c)')
                
        if chi[c+1].type == 1:
            if chi[c+1].b1 > chi[c+1].b-2*dat.b1_ext_in[2] or chi[c+1].c_t > chi[c+1].c_m or chi[c+1].c_m > chi[c+1].c_r:
                raise TypeError('Problema na geometria da planta (c+1)')
        
        # Zerar a pontuação dos filhos 
        chi[c].score = 0
        chi[c+1].score = 0
        
        c += 2
    
    
    # Mutação
    select1 = np.array([int(x) for x in np.random.uniform(0,1,len(pop)) <= dat.mu])
    if sum(select1) != 0:
        print('<< Mutação >>')
        select2 = np.argwhere(select1 == 1)
        k = 0
        for i in select2:
            s = select2[k,0]
            
            
            # A mutação funciona de modo a atribuir novos valores aos respectivos campos
            
            # Genes a serem alterados 
            # Envergadura b
            # Envergadura b1 (asas bitrapezoidais apenas)
            # Corda da raiz
            # Corda do meio (asas bitrapezoidais apenas)
            # Corda da ponta
            # - aerofólio da raiz
            # - aerofólio do meio (asas bitrapezoidais apenas)
            # - aerofólio da ponta
            # Torção geométrica na ponta
            
            mu_op_v = np.random.randint(2,size=12)
            
            # Envergadura b
            if mu_op_v[0] == 1:
                chi[s].b = dat.b_ext[random.randint(0,len(dat.b_ext)-1)]
            
            # Envergadura da primeira seção b1 (asas bitrapezoidais apenas)
            if mu_op_v[1] == 1 and chi[s].type == 1:
                if dat.b1_ext_in[1] >= chi[s].b - 2*dat.b1_ext_in[2]:
                    dat.b1_ext = np.arange(dat.b1_ext_in[0],chi[s].b-dat.b1_ext_in[2]*2+dat.b1_step,dat.b1_step)
                else:
                    dat.b1_ext = np.arange(dat.b1_ext_in[0],dat.b1_ext_in[1]+dat.b1_step,dat.b1_step)
                
                chi[s].b1 = dat.b1_ext[random.randint(0,len(dat.b1_ext)-1)]
            
            # Corda da raiz c_r
            if mu_op_v[2] == 1:
                chi[s].c_r = dat.c_r_ext[random.randint(0,len(dat.c_r_ext)-1)]
                
                # Se a corda da ponta for maior do que a corda da raiz, atribuir
                # o valor da ponta à raiz
                if chi[s].c_t > chi[s].c_r:
                    chi[s].c_r = chi[s].c_t
            
            # Corda da ponta c_t (deve cumprir o requisito c_t <= c_r)
            if mu_op_v[3] == 1:
                if chi[s].c_r < dat.c_t_ext_in[1]:
                    dat.c_t_ext = np.arange(dat.c_t_ext_in[0],chi[s].c_r+dat.c_t_step,dat.c_t_step)
                else:
                    dat.c_t_ext = np.arange(dat.c_t_ext_in[0],dat.c_t_ext_in[1]+dat.c_t_step,dat.c_t_step)
                
                chi[s].c_t = dat.c_t_ext[random.randint(0,len(dat.c_t_ext)-1)]     
                
                # Se a corda da ponta for maior do que a corda da raiz, atribuir
                # o valor da raiz à ponta
                if chi[s].c_t > chi[s].c_r:
                    chi[s].c_t = chi[s].c_r
            
            # Corda do meio c_m (asas bitrapezoidais apenas)
            if mu_op_v[4] == 1 and not 'L' in dat.c_m_ext_in:
                dat.c_m_ext = [0,0]
                if dat.c_m_ext_in[0] < chi[s].c_t:
                    # Caso o valor mínimo da extensão da corda do meio seja menor do que
                    # a corda da ponta da asa atual, usar o valor da corda da ponta da
                    # asa como o valor mínimo da extensão
                    dat.c_m_ext[0] = chi[s].c_t
                else:
                    # Caso contrário, usar o valor original
                    dat.c_m_ext[0] = dat.c_m_ext_in[0]
                
                if dat.c_m_ext_in[1] > chi[s].c_t:
                    # Caso o valor máximo da extensão da corda do meio seja maior do que
                    # a corda da raiz da asa atual, usar o valor da corda da raiz da
                    # asa como o valor mínimo da extensão
                    dat.c_m_ext[1] = chi[s].c_r
                else:
                    # Caso contrário, usar o valor original
                    dat.c_m_ext[1] = dat.c_m_ext_in[1]
                
                chi[s].c_m = dat.c_m_ext[random.randint(0,len(dat.c_m_ext)-1)]
            
            # Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi[s].type == 1:
                
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[s].b1 > chi[s].b-dat.b1_ext_in[2]*2:
                    chi[s].b1 = chi[s].b-dat.b1_ext_in[2]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[s].c_m > chi[s].c_r and not 'L' in dat.c_m_ext_in:
                    chi[s].c_m = chi[s].c_r
                
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[s].c_m < chi[s].c_t and not 'L' in dat.c_m_ext_in:
                    chi[s].c_m = chi[s].c_t

            # Aerofólio da raiz
            if mu_op_v[5] == 1:
                
                op = np.random.randint(2,size=3)
                
                if op[0] == 1:
                    chi[s].af_r[0] = random.randint(dat.m_ext_r[0],dat.m_ext_r[1])
                    # Terminar o resto da transformação para perfil simétrico
                    if chi[s].af_r[0] == 0:
                        chi[s].af_r[1] = 0
                
                if ( op[1] == 1 and chi[s].af_r[0] != 0 ) and ( chi[s].af_r[0] != 0 and chi[s].af_r[1] == 0 ):
                    # Considerar também o resto da transformação simétrico -> assimétrico
                    chi[s].af_r[1] =  random.randint(dat.p_ext_r[0],dat.p_ext_r[1])
                
                if op[2] == 1:
                    chi[s].af_r[2] = random.randint(dat.t_ext_r[0],dat.t_ext_r[1])
            
            # Aerofólio do meio (asas bitrapezoidais apenas)
            if mu_op_v[6] == 1 and chi[s].type == 1:
                op = np.random.randint(2,size=3)
                
                if op[0] == 1:
                    chi[s].af_m[0] = random.randint(dat.m_ext_m[0],dat.m_ext_m[1])
                    # Terminar o resto da transformação para perfil simétrico
                    if chi[s].af_m[0] == 0:
                        chi[s].af_m[1] = 0
                
                if ( op[1] == 1 and chi[s].af_m[0] != 0 ) and ( chi[s].af_m[0] != 0 and chi[s].af_m[1] == 0 ):
                    # Considerar também o resto da transformação simétrico -> assimétrico
                    chi[s].af_m[1] =  random.randint(dat.p_ext_m[0],dat.p_ext_m[1])
                
                if op[2] == 1:
                    chi[s].af_m[2] = random.randint(dat.t_ext_m[0],dat.t_ext_m[1])
                    
            # Aerofólio da ponta
            if mu_op_v[7] == 1:
                op = np.random.randint(2,size=3)
                
                if op[0] == 1:
                    chi[s].af_t[0] = random.randint(dat.m_ext_t[0],dat.m_ext_t[1])
                    # Terminar o resto da transformação para perfil simétrico
                    if chi[s].af_t[0] == 0:
                        chi[s].af_t[1] = 0
                
                if ( op[1] == 1 and chi[s].af_t[0] != 0 ) and ( chi[s].af_t[0] != 0 and chi[s].af_t[1] == 0 ):
                    # Considerar também o resto da transformação simétrico -> assimétrico
                    chi[s].af_t[1] =  random.randint(dat.p_ext_t[0],dat.p_ext_t[1])
                
                if op[2] == 1:
                    chi[s].af_t[2] = random.randint(dat.t_ext_t[0],dat.t_ext_t[1])
            
            # Torção geométrica na ponta
            if mu_op_v[8] == 1:
                chi[s].tw_t = dat.tw_t_ext[random.randint(0,len(dat.tw_t_ext)-1)]
            
            # Enflechamento total (asas trapezoidais simples)
            if mu_op_v[9] == 1 and chi[s].type == 0 and not 'Z' in dat.sweep_ext_in:
                chi[s].sweep = dat.sweep_ext[random.randint(0,len(dat.sweep_ext)-1)]
            
            # Enflechamento da primeira seção (asas trapezoidais duplas)
            if mu_op_v[10] == 1 and chi[s].type == 1 and not 'Z' in dat.sweep1_ext_in:
                chi[s].sweep1 = dat.sweep1_ext[random.randint(0,len(dat.sweep1_ext)-1)]
            
            # Enflechamento da segunda seção (asas trapezoidais duplas)
            if mu_op_v[11] == 1 and chi[s].type == 1 and not 'Z' in dat.sweep2_ext_in:
                chi[s].sweep2 = dat.sweep2_ext[random.randint(0,len(dat.sweep2_ext)-1)]
            
            # Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
            if chi[s].type == 1:
                if chi[s].b1 > chi[s].b-2*dat.b1_ext_in[2] or chi[s].c_t > chi[s].c_m and not 'L' in dat.c_m_ext_in or chi[s].c_m > chi[s].c_r and not 'L' in dat.c_m_ext_in:
                    raise TypeError('Problema na geometria da planta')
            
            k += 1

    # substituir a população inicial pelos filhos
    for i in range(dat.N):
        pop[i] = chi[i].deepcopy()
    
	# Aplicar elitismo
    if dat.elite == 1:
        # Passar o melhor indivíduo pra nova população
        pop[0] = archive[loop].deepcopy()
        pop[0].score = 0


# Final ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end = time()
mi = np.fix((end-start)/60)
s = (end-start)%60
print('Tempo: ' + str(mi) + 'min e ' + str(s) + 's')

# Fazer gráficos dos coeficientes dos melhores indivíduos de cada geração
for i in range(dat.cases):
    aero_m = np.zeros((dat.iter,4))
    for j in range(dat.iter):
        aero_m[j,:] = archive[j].aero[i,:]
    plt.figure()
    plt.subplots()
    plt.plot(np.arange(1,dat.iter+1),aero_m[:,0],'g-*',label='CL')
    plt.plot(np.arange(1,dat.iter+1),aero_m[:,1],'r-*',label='CD')
    plt.plot(np.arange(1,dat.iter+1),aero_m[:,3],'b-*',label='CM')
    plt.ylabel('CL, CD, CM')
    plt.legend(loc='upper left')
    plt.twinx()
    plt.plot(np.arange(1,dat.iter+1),aero_m[:,2],'k-*',label='L/D')
    plt.ylabel('L/D')
    plt.legend(loc='upper right')
    plt.grid('True')
    plt.title('Melhores resultados - Condição de voo '  + str(i+1) + ':  Re ' + str(dat.reynolds[i]) + ', AoA ' + str(dat.aoa[i]) + '°')
    plt.xlabel('Iteração')  
    
## Pegar o struct de arquivo e imprimir todos

for i in range(len(archive)):
    print('<< Iteração {} >>'.format(i+1))
    print(' - Dados da planta - ')
    if archive[i].type == 0:
        print('Tipo: trapezoidal simples (type = 0)')
    else:
        print('Tipo: trapezoidal dupla (type = 1)')
    print('Envergadura b = ' + str(archive[i].b))
    if archive[i].type == 1:print('Envergadura da primeira seção b1 = ' + str(archive[i].b1))
    print('Corda da raiz c_r = ' + str(archive[i].c_r))
    if archive[i].type == 1:print('Corda do meio c_m = ' + str(archive[i].c_m))
    print('Corda da ponta c_t = ' + str(archive[i].c_t))
    if archive[i].type == 0:
        print('Enflechamento completo sweep = ' + str(archive[i].sweep) + '°')
    else:
        print('Enflechamento da primeira seção sweep1 = ' + str(archive[i].sweep1) + '°')
        print('Enflechamento da segunda seção sweep2 = ' + str(archive[i].sweep2) + '°')
    print('')
    
    
    print(' - Dados dos aerofólios - ')
    
    print('Raiz: NACA ' + str(int(archive[i].af_r[0])) + str(int(archive[i].af_r[1])) + str(int(archive[i].af_r[2])))
    if archive[i].type == 1: print('Meio: NACA ' + str(int(archive[i].af_m[0])) + str(int(archive[i].af_m[1])) + str(int(archive[i].af_m[2])))
    print('Ponta: NACA ' + str(int(archive[i].af_t[0])) + str(int(archive[i].af_t[1])) + str(int(archive[i].af_t[2])))
    print('Torção geométrica na ponta tw_t = ' + str(archive[i].tw_t))
    
    
    print(' - Dados aerodinâmicos - ')
    for j in range(dat.cases):
        print('Condição de voo %d: %f, %f, %f, %f'%(j,archive[i].aero[j,0],archive[i].aero[j,1],archive[i].aero[j,2],archive[i].aero[j,3]))
        
        if 'q' in dat.coeff_op[:,0] or '#' in dat.coeff_op[:,0]:
            print('L = %f N'%(archive[i].aero[j,0]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))
        
        if 'q' in dat.coeff_op[:,1] or '#' in dat.coeff_op[:,1]:
            print('D = %f N'%(archive[i].aero[j,1]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))
        
        if 'q' in dat.coeff_op[:,3]:
            print('M = %f Nm'%(archive[i].aero[j,3]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))

    print('Pontuação: %f'%(archive[i].score))
    print('')
    print('')
    