## Algoritmo genético
# Otimização de asas trapezoidais simples ou duplas, com perfis CST, feita
# em termos de uma geometria inicial específica

# Pacotes necessários
import numpy as np
from ypstruct import structure
import random
# import os
import matplotlib.pyplot as plt
from time import time
# Funções utilizadas
from error_check_cst_3D_SPECIFY import error_check_cst_3D_SPECIFY
from quality import quality
from run_cst_TCC2_3D import run_cst_TCC2_3D
from run_apame_mesher_cst import run_apame_mesher_cst
from run_apame import run_apame
# from run_openvsp_cst import run_openvsp_cst
from run_openvsp_cst_VLM_specify import run_openvsp_cst_VLM_specify
# from run_vspaero import run_vspaero
from run_vspaero_VLM import run_vspaero_VLM
from fitness_cst_3D import fitness_cst_3D
# from make_vector import make_vector
from selection_crossover import selection_crossover
from plot_planform import plot_planform 
from plot_airfoil_cst_TCC2 import plot_airfoil_cst_TCC2

start = time()

# Parâmetros do algoritmo
dat = structure()
dat.N = 20                              # Número de indivíduos na população
dat.mu = 0.05                           # Probabilidade de mutação (definida entre zero e um)
dat.iter = 5                           # Número de iterações
dat.elite = 1                         # Aplicar elitismo?
dat.subs = 1                          # Substituir asas sem resultados? (ver ainda se isto será necessário)

# Dados da asa a ser otimizada
# Planta
dat.type = 0 # (0 -> trapezoidal simples, 1 -> trapezoidal dupla) 
dat.or_b = 4.563
dat.or_b1 = 10
dat.or_c_r = 2
dat.or_c_m = 1
dat.or_c_t = 0.34
dat.or_sweep = 35
dat.or_sweep1 = 'Z'
dat.or_sweep2 = 'Z'
# Aerofólio da raiz (NACA 65-006)
dat.or_v_ex_r = [0.001964, 0.062026, 0.101569, 0.088828, 2.086990, 0.000000]
dat.or_v_in_r = [0.001964, 0.062026, 0.101569, 0.088828, 2.086990, 0.000000]
dat.symm_override_r = 0
# Aerofólio do meio
dat.or_v_ex_m = [0.01, 0.2, 0.2, 0.1, 10, 0]
dat.or_v_in_m = [0.01, 0.1, 0.1, 0.1, 10, 0]
dat.symm_override_m = 0;
# Aerofólio da ponta (NACA 65-006)
dat.or_v_ex_t = [0.001964, 0.062026, 0.101569, 0.088828, 2.086990, 0.000000]
dat.or_v_in_t = [0.001964, 0.062026, 0.101569, 0.088828, 2.086990, 0.000000]
dat.symm_override_t = 0
dat.or_tw_t = 0

# Parâmetros da geometria: planta da asa
#dat.planf_op = 0.5; # Proporção de asas trapezoidais simples e bitrapezoidais (0->todas trapezoidais simples, 1->todas bitrapezoidais)
dat.b_ext = [-1,1,5] # Envergadura completa [m] (limite inferior, valor mínimo possível)
dat.b1_ext = [-0.2,0.2,1,0.5] # Envergadura da raiz ao meio [m] (asas bitrapezoidais apenas) (valor mínimo,valor máximo,valor mínimo possível,separação mínima da ponta da asa (considerando apenas uma metade)) 
dat.c_r_ext = [-0.5,0.5,0.5] # Corda da raiz [m] (limite inferior, valor mínimo possível)
dat.c_m_ext = [-0.2,0.2,0.5] # Corda do meio [m] (asas bitrapezoidais apenas) (limite inferior, limite superior, valor mínimo possível) (a opção 'L' força o formato trapezoidal simples)
dat.c_t_ext = [-0.5,0.5,0.2] # Corda da ponta [m] (limite inferior, valor mínimo possível)
dat.sweep_ext = [-0,0] # Enflechamento de asas trapezoidais simples (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep1_ext = [-0,0] # Enflechamento da primeira seção de asas trapezoidais duplas (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep2_ext = [-0,0] # Enflechamento da segunda seção de asas trapezoidais duplas (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
# Parâmetros da geometria: aerofólio da raiz
dat.BPn_r = len(dat.or_v_ex_r)-2                  # Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_r = 0.5
dat.N2_r = 1
dat.le_R_ext1_r = [-0.01,0.01,0.005] # Limite inferior, limite superior, valor mínimo possível (pois Rle > 0)
dat.le_R_ext2_r = [-0.01,0.01,0.005]
dat.A_ext1_r = [-0.1,0.1]
dat.A_ext2_r = [-0.1,0.1]
dat.B_ext1_r = [-5,5] # Limites inferior e superior
dat.B_ext2_r = [dat.or_v_ex_r[-2] + dat.or_v_in_r[-2],5] # O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
# Parâmetros da geometria: aerofólio do meio (asas bitrapezoidais apenas) 
dat.BPn_m = len(dat.or_v_ex_m)-2                  # Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_m = 0.5
dat.N2_m = 1
dat.le_R_ext1_m = [-0.01,0.01,0.005] # (admite opção 'L')
dat.le_R_ext2_m = [-0.01,0.01,0.005]
dat.A_ext1_m = [-0.1,0.1]
dat.A_ext2_m = [-0.1,0.1]
dat.B_ext1_m = [-5,5] # Limites inferior e superior
dat.B_ext2_m = [dat.or_v_ex_m[-2] + dat.or_v_in_m[-2],5] # O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
# Parâmetros da geometria: aerofólio da ponta
dat.BPn_t = len(dat.or_v_ex_t)-2                  # Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_t = 0.5
dat.N2_t = 1
dat.le_R_ext1_t = [-0.01,0.01,0.005]
dat.le_R_ext2_t = [-0.01,0.01,0.005]
dat.A_ext1_t = [-0.1,0.1]
dat.A_ext2_t = [-0.1,0.1]
dat.B_ext1_t = [-5,-5] # Limites inferior e superior
dat.B_ext2_t = [dat.or_v_ex_t[-2] + dat.or_v_in_t[-2],5] # O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
dat.tw_t_ext = [-1,0]

# Torção geométrica
dat.tw_t_ext_in = [-1,1] # Torção geométrica [°]

# Parâmetros da malha
dat.np = 30 # Número de pontos na geração de ordenadas nos aerofólios
dat.np_op = 1 # 1 -> cosspace, 0 -> cosspace_half
dat.nb = [1,1] # Número de seções intermediárias (raiz/ponta) [número de seções,0] ou [concentração por metro,1]
dat.nb1 = 'L' # Número de seções intermediárias (raiz/meio) (asas bitrapezoidais apenas) (opção 'L' faz com que nb1 e nb2 sejam uniformemente determinados ao longo da envergadura)
dat.nb2 = None # Número de seções intermediárias (meio/ponta) (asas bitrapezoidais apenas)

# Parâmetros das simulações
dat.method = 2                        # Método da simulação (1->painéis 3D no APAME, 2->VLM no VSPAERO)
dat.cases = 1                          # Número de condições de voo a serem analisadas
dat.v_ref = [100,100,100] # Velocidades de referência [m/s] 
dat.rho = [1.225,1.225,1.225] # Densidades do ar [kg/m^3] 
dat.p_atm = [101325,101325,101325] # Pressões do ar [Pa] (irrelevante neste algoritmo)
dat.mach = [0.6,0,0.] # Números de Mach
dat.reynolds = [1e6,1e6,1e6]           # Valores dos números de Reynolds para as simulações (irrelevante neste algoritmo))
dat.aoa = [2,0,4]                      # Ângulos de ataque
dat.karman_tsien = ['Y','N','N']    # Correção de compressibilidade de Karman-Tsien
dat.wake_iters = [1,1,1]            # Números de iterações no VSPAERO
dat.coeff_op = np.array((['!','^','^','!'],       # Uma linha para cada condição de voo
                         ['!','!','!','!'],
                         ['!','!','!','!']))
dat.coeff_val = np.array(([0.05,7e-3,90,-1e-1],
                 [0.5,0,0,-0.08],
                 [0,0,0,-0.08]))
dat.coeff_F = np.array(([1,1,1,1],
               [1,1,1,1],
               [1,1,1,1]))
# [CL CD L/D CM] Definição de cada linha da matriz dat.coeff_op
# '!' -> não usar como função objetivo
# '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD)
# 'c' -> buscar valor constante de coeficiente de momento (arbitrário)
# 'k' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val(1,4))
# 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
# em sua respectiva casa na matriz dat.coeff_val
# 'q' -> procurar por um valor específico de força de sustentação, força de arrasto
# ou momento de arfagem (CL, CD e CM). Nesse caso, definir o valor em sua casa na matriz dat.coeff_val
# '#' -> procurar por um valor máximo (L) ou mínimo (D)
# A matriz dat.coeff_F dá os pesos de cada função objetivo

# Checagem de erros
dat = error_check_cst_3D_SPECIFY(dat)

# Template dos structs
empty = structure()
# Forma da planta
empty.type = dat.type
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
# Aerofólios 
empty.v_ex_r = np.zeros(dat.BPn_r+2)
empty.v_in_r = np.zeros(dat.BPn_r+2)
empty.symm_r = None
empty.v_ex_m = np.zeros(dat.BPn_m+2)
empty.v_in_m = np.zeros(dat.BPn_m+2)
empty.symm_m = None
empty.v_ex_t = np.zeros(dat.BPn_t+2)
empty.v_in_t = np.zeros(dat.BPn_t+2)
empty.symm_t = None
empty.tw_m = None
empty.tw_t = None
# Dados da malha
empty.NODE = None
empty.PANEL = None
empty.sec_N = None
# Dados aerodinâmicos e pontuação
empty.aero = None # Terá o mesmo formato que a matriz coeff_op
empty.score = 0

# Inicializar os structs
pop = empty.repeat(dat.N)
chi = empty.repeat(dat.N)

# Fazer um struct da asa original
original = empty.deepcopy()
original.type = dat.type
original.b = dat.or_b
original.b1 = dat.or_b1
original.c_r = dat.or_c_r
original.c_m = dat.or_c_m
original.c_t = dat.or_c_t
original.v_ex_r = dat.or_v_ex_r
original.sweep = dat.or_sweep
original.sweep1 = dat.or_sweep1
original.sweep2 = dat.or_sweep2
if dat.symm_override_r == 1:
    original.v_in_r = dat.or_v_ex_r
else:
    original.v_in_r = dat.or_v_in_r
original.v_ex_m = dat.or_v_ex_m
if dat.symm_override_m == 1:
    original.v_in_m = dat.or_v_ex_m
else:
    original.v_in_m = dat.or_v_in_m
original.v_ex_t = dat.or_v_ex_t
if dat.symm_override_t == 1:
    original.v_in_t = dat.or_v_ex_t
else:
    original.v_in_t = dat.or_v_in_t
original.tw_m = 'L'
original.tw_t = dat.or_tw_t
# Pegar dados aerodinâmicos da asa original para fins de comparação
if dat.method == 1: # Painéis 3D
    original = run_apame_mesher_cst(original,dat,1) # Obter malha do apame
    original.aero = run_apame(original,dat) # Fazer simulação
elif dat.method == 2: # VLM
    original = run_openvsp_cst_VLM_specify(original,dat) # Obter geometria OpenVSP            
    original.aero = run_vspaero_VLM(original,dat) # Simulação
    
if original.aero == 'n':
    i = input('A simulação do indivíduo original não convergiu. Continuar mesmo assim? (y/n) ')
    if i != 'y' and i != 'Y':
        raise TypeError('Execução terminada')

# Gerar população inicial
for i in range(dat.N):
    print('Indivíduo ' + str(i+1))
    
    # Gerar forma da planta
	# Envergadura
    pop[i].b = dat.or_b + random.random()*dat.b_ext[random.randint(0,1)]
    if pop[i].b < dat.b_ext[2]:
        pop[i].b = dat.b_ext[2]
    
	# Corda da raiz
    pop[i].c_r = dat.or_c_r + random.random()*dat.c_r_ext[random.randint(0,1)]
    if pop[i].c_r < dat.c_r_ext[2]:
        pop[i].c_r = dat.c_r_ext[2]

	# Corda da ponta
    pop[i].c_t = dat.or_c_t + random.random()*dat.c_t_ext[random.randint(0,1)]
    if pop[i].c_t < dat.c_t_ext[2]: # Se for menor que o mínimo permitido
        pop[i].c_t = dat.c_t_ext[2]
    elif pop[i].c_t > pop[i].c_r: # Se for maior que a corda da raiz
        pop[i].c_t = pop[i].c_r
    
	# Torção geométrica na ponta
    pop[i].tw_t = dat.or_tw_t + random.random()*dat.tw_t_ext[random.randint(0,1)]
    
    # Gerar aerofólio da raiz
    check = 0
    while check == 0:
        
        if dat.or_v_ex_r == dat.or_v_in_r or dat.symm_override_r == 1: # Perfil simétrico
            # Vetor com informações do extradorso
            pop[i].v_ex_r[0] = dat.or_v_ex_r[0] + random.random()*dat.le_R_ext1_r[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn_r):
                pop[i].v_ex_r[a] = dat.or_v_ex_r[a] + random.random()*dat.A_ext1_r[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_ex_r[dat.BPn_r] = dat.or_v_ex_r[dat.BPn_r] + random.random()*dat.B_ext1_r[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex_r[dat.BPn_r+1] = dat.or_v_ex_r[dat.BPn_r+1] #randi(dat.delta_range)*0.1*rand   # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in_r = pop[i].v_ex_r
        
            # Checar o raio do bordo de ataque
            if pop[i].v_ex_r[0] < dat.le_R_ext1_r[2]:
                pop[i].v_ex_r[0] = dat.le_R_ext1_r[2]
                pop[i].v_in_r[0] = pop[i].v_ex_r[0]
            
            # Checar separação do bordo de fuga
            if pop[i].v_ex_r[dat.BPn_r] + pop[i].v_in_r[dat.BPn_r] < dat.B_ext2_r[0]:
                pop[i].v_ex_r[dat.BPn_r] = dat.B_ext2_r[0]/2
                pop[i].v_in_r[dat.BPn_r] = dat.B_ext2_r[0]/2
            
            pop[i].symm_r = 1
            
        else: # perfil assimétrico
            # Vetor com informações do extradorso
            pop[i].v_ex_r[0] = dat.or_v_ex_r[0] + random.random()*dat.le_R_ext1_r[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn_r):
                pop[i].v_ex_r[a] = dat.or_v_ex_r[a] + random.random()*dat.A_ext1_r[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_ex_r[dat.BPn_r] = dat.or_v_ex_r[dat.BPn_r] + random.random()*dat.B_ext1_r[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex_r[dat.BPn_r+1] = dat.or_v_ex_r[dat.BPn_r+1]  # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in_r[0] = dat.or_v_in_r[0]+ random.random()*dat.le_R_ext2_r[random.randint(0,1)]                     # Raio do bordo de ataque
            for a in range(1,dat.BPn_r):
                pop[i].v_in_r[a] = dat.or_v_in_r[a] + random.random()*dat.A_ext2_r[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_in_r[dat.BPn_r] = (dat.B_ext2_r[0] - pop[i].v_ex_r[dat.BPn_r]) + random.random()*dat.B_ext2_r[1]     # Ângulo do bordo de fuga
            pop[i].v_in_r[dat.BPn_r+1] = dat.or_v_in_r[dat.BPn_r+1] 
            
            # Checar o raio do bordo de ataque
            if pop[i].v_ex_r[0] < dat.le_R_ext1_r[2]:
                pop[i].v_ex_r[0] = dat.le_R_ext1_r[2]
        
            if pop[i].v_in_r[0] < dat.le_R_ext2_r[2]:
                pop[i].v_in_r[0] = dat.le_R_ext2_r[2]
            
            # Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop[i].v_in_r[dat.BPn_r] < dat.B_ext2_r[0] - pop[i].v_ex_r[dat.BPn_r]:
                pop[i].v_in_r[dat.BPn_r] = dat.B_ext2_r[0] - pop[i].v_ex_r[dat.BPn_r]
            
            # Checar os pesos (soma de pesos do intradorso deve ser menor ou
            # igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop[i].v_ex_r[1:dat.BPn_r])
            sum2 = sum(pop[i].v_in_r[1:dat.BPn_r])
            if sum2 > sum1:continue
            
            pop[i].symm_r = 0
        
        # Checagem de qualidade
        check = quality(run_cst_TCC2_3D(pop[i].v_ex_r,pop[i].v_in_r,dat,[dat.N1_r,dat.N2_r]),dat)
        
    
    # Gerar aerofólio da ponta
    check = 0
    while check == 0:
        
        if dat.or_v_ex_t == dat.or_v_in_t or dat.symm_override_t == 1: # Perfil simétrico
            # Vetor com informações do extradorso
            pop[i].v_ex_t[0] = dat.or_v_ex_t[0] + random.random()*dat.le_R_ext1_t[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn_t):
                pop[i].v_ex_t[a] = dat.or_v_ex_t[a] + random.random()*dat.A_ext1_t[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_ex_t[dat.BPn_t] = dat.or_v_ex_t[dat.BPn_t] + random.random()*dat.B_ext1_t[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex_t[dat.BPn_t+1] = dat.or_v_ex_t[dat.BPn_t+1] #randi(dat.delta_tange)*0.1*rand   # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in_t = pop[i].v_ex_t
        
            # Checar o raio do bordo de ataque
            if pop[i].v_ex_t[0] < dat.le_R_ext1_t[2]:
                pop[i].v_ex_t[0] = dat.le_R_ext1_t[2]
                pop[i].v_in_t[0] = pop[i].v_ex_t[0]
            
            # Checar separação do bordo de fuga
            if pop[i].v_ex_t[dat.BPn_t] + pop[i].v_in_t[dat.BPn_t] < dat.B_ext2_t[0]:
                pop[i].v_ex_t[dat.BPn_t] = dat.B_ext2_t[0]/2
                pop[i].v_in_t[dat.BPn_t] = dat.B_ext2_t[0]/2
            
            pop[i].symm_t = 1
            
        else: # perfil assimétrico
            # Vetor com informações do extradorso
            pop[i].v_ex_t[0] = dat.or_v_ex_t[0] + random.random()*dat.le_R_ext1_t[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn_t):
                pop[i].v_ex_t[a] = dat.or_v_ex_t[a] + random.random()*dat.A_ext1_t[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_ex_t[dat.BPn_t] = dat.or_v_ex_t[dat.BPn_t] + random.random()*dat.B_ext1_t[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex_t[dat.BPn_t+1] = dat.or_v_ex_t[dat.BPn_t+1]  # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in_t[0] = dat.or_v_in_t[0]+ random.random()*dat.le_R_ext2_t[random.randint(0,1)]                     # Raio do bordo de ataque
            for a in range(1,dat.BPn_t):
                pop[i].v_in_t[a] = dat.or_v_in_t[a] + random.random()*dat.A_ext2_t[random.randint(0,1)]        # Pesos intermediários
            
            pop[i].v_in_t[dat.BPn_t] = (dat.B_ext2_t[0] - pop[i].v_ex_t[dat.BPn_t]) + random.random()*dat.B_ext2_t[1]     # Ângulo do bordo de fuga
            pop[i].v_in_t[dat.BPn_t+1] = dat.or_v_in_t[dat.BPn_t+1] 
            
            # Checar o raio do bordo de ataque
            if pop[i].v_ex_t[0] < dat.le_R_ext1_t[2]:
                pop[i].v_ex_t[0] = dat.le_R_ext1_t[2]
        
            if pop[i].v_in_t[0] < dat.le_R_ext2_t[2]:
                pop[i].v_in_t[0] = dat.le_R_ext2_t[2]
            
            # Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop[i].v_in_t[dat.BPn_t] < dat.B_ext2_t[0] - pop[i].v_ex_t[dat.BPn_t]:
                pop[i].v_in_t[dat.BPn_t] = dat.B_ext2_t[0] - pop[i].v_ex_t[dat.BPn_t]
            
            # Checar os pesos (soma de pesos do intradorso deve ser menor ou
            # igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop[i].v_ex_t[1:dat.BPn_t])
            sum2 = sum(pop[i].v_in_t[1:dat.BPn_t])
            if sum2 > sum1:continue
            
            pop[i].symm_t = 0
        
        # Checagem de qualidade
        check = quality(run_cst_TCC2_3D(pop[i].v_ex_t,pop[i].v_in_t,dat,[dat.N1_t,dat.N2_t]),dat)
    
    # Dados adicionais para asas bitrapezoidais
    if pop[i].type == 1:
        # Mais dados da planta
        if dat.or_c_m == 'L':
            # Caso c_m seja definido como 'L', seu valor real será atribúido pela
            # função run_apame_mesher_cst
            pop[i].c_m = 'L'
        else:
            # Caso contrário, é realizado o processo abaixo
            pop[i].c_m = dat.or_c_m + random.random()*dat.c_m_ext[random.randint(0,1)]
            if pop[i].c_m < dat.c_m_ext[2]:
                pop[i].c_m = dat.c_m_ext[2]
                
        pop[i].b1 = dat.or_b1 + random.random()*dat.b1_ext[random.randint(0,1)]
        if pop[i].b1 < dat.b1_ext[2]:
            pop[i].b1 = dat.b1_ext[2]
        
        # Aerofólio do meio
        pop[i].tw_m = 'L' # Essa sempre será a configuração deste algoritmo, mas isso pode ser alterado (com as devidas alterações no resto do código)
        check = 0
        while check == 0:
            
            if dat.or_v_ex_m == dat.or_v_in_m or dat.symm_override_m == 1: # Perfil simétrico
                # Vetor com informações do extradorso
                pop[i].v_ex_m[0] = dat.or_v_ex_m[0] + random.random()*dat.le_R_ext1_m[random.randint(0,1)]    # Raio do bordo de ataque
                for a in range(1,dat.BPn_m):
                    pop[i].v_ex_m[a] = dat.or_v_ex_m[a] + random.random()*dat.A_ext1_m[random.randint(0,1)]        # Pesos intermediários
                
                pop[i].v_ex_m[dat.BPn_m] = dat.or_v_ex_m[dat.BPn_m] + random.random()*dat.B_ext1_m[random.randint(0,1)]     # Ângulo do bordo de fuga
                pop[i].v_ex_m[dat.BPn_m+1] = dat.or_v_ex_m[dat.BPn_m+1] #randi(dat.delta_mange)*0.1*rand   # delta_z
                
                # Vetor com informações do intradorso
                pop[i].v_in_m = pop[i].v_ex_m
            
                # Checar o raio do bordo de ataque
                if pop[i].v_ex_m[0] < dat.le_R_ext1_m[2]:
                    pop[i].v_ex_m[0] = dat.le_R_ext1_m[2]
                    pop[i].v_in_m[0] = pop[i].v_ex_m[0]
                
                # Checar separação do bordo de fuga
                if pop[i].v_ex_m[dat.BPn_m] + pop[i].v_in_m[dat.BPn_m] < dat.B_ext2_m[0]:
                    pop[i].v_ex_m[dat.BPn_m] = dat.B_ext2_m[0]/2
                    pop[i].v_in_m[dat.BPn_m] = dat.B_ext2_m[0]/2
                
                pop[i].symm_m = 1
                
            else: # perfil assimétrico
                # Vetor com informações do extradorso
                pop[i].v_ex_m[0] = dat.or_v_ex_m[0] + random.random()*dat.le_R_ext1_m[random.randint(0,1)]    # Raio do bordo de ataque
                for a in range(1,dat.BPn_m):
                    pop[i].v_ex_m[a] = dat.or_v_ex_m[a] + random.random()*dat.A_ext1_m[random.randint(0,1)]        # Pesos intermediários
                
                pop[i].v_ex_m[dat.BPn_m] = dat.or_v_ex_m[dat.BPn_m] + random.random()*dat.B_ext1_m[random.randint(0,1)]     # Ângulo do bordo de fuga
                pop[i].v_ex_m[dat.BPn_m+1] = dat.or_v_ex_m[dat.BPn_m+1]  # delta_z
                
                # Vetor com informações do intradorso
                pop[i].v_in_m[0] = dat.or_v_in_m[0]+ random.random()*dat.le_R_ext2_m[random.randint(0,1)]                     # Raio do bordo de ataque
                for a in range(1,dat.BPn_m):
                    pop[i].v_in_m[a] = dat.or_v_in_m[a] + random.random()*dat.A_ext2_m[random.randint(0,1)]        # Pesos intermediários
                
                pop[i].v_in_m[dat.BPn_m] = (dat.B_ext2_m[0] - pop[i].v_ex_m[dat.BPn_m]) + random.random()*dat.B_ext2_m[1]     # Ângulo do bordo de fuga
                pop[i].v_in_m[dat.BPn_m+1] = dat.or_v_in_m[dat.BPn_m+1] 
                
                # Checar o raio do bordo de ataque
                if pop[i].v_ex_m[0] < dat.le_R_ext1_m[2]:
                    pop[i].v_ex_m[0] = dat.le_R_ext1_m[2]
            
                if pop[i].v_in_m[0] < dat.le_R_ext2_m[2]:
                    pop[i].v_in_m[0] = dat.le_R_ext2_m[2]
                
                # Checar separação do bordo de fuga (beta2>=L-beta1)
                if pop[i].v_in_m[dat.BPn_m] < dat.B_ext2_m[0] - pop[i].v_ex_m[dat.BPn_m]:
                    pop[i].v_in_m[dat.BPn_m] = dat.B_ext2_m[0] - pop[i].v_ex_m[dat.BPn_m]
                
                # Checar os pesos (soma de pesos do intradorso deve ser menor ou
                # igual à soma de pesos do extradorso - pesos intermediários)
                sum1 = sum(pop[i].v_ex_m[1:dat.BPn_m])
                sum2 = sum(pop[i].v_in_m[1:dat.BPn_m])
                if sum2 > sum1:continue
                
                pop[i].symm_m = 0
            
            # Checagem de qualidade
            check = quality(run_cst_TCC2_3D(pop[i].v_ex_m,pop[i].v_in_m,dat,[dat.N1_m,dat.N2_m]),dat)
       
        # Aplicar checagens de geometria da planta
        if pop[i].c_m > pop[i].c_r and dat.or_c_m != 'L':
            pop[i].c_m = pop[i].c_r
        
        if pop[i].c_t > pop[i].c_m and dat.or_c_m != 'L':
            pop[i].c_t = pop[i].c_m
        
        if pop[i].b1 > pop[i].b - 2*dat.b1_ext[3]:
            pop[i].b1 = pop[i].b - 2*dat.b1_ext[3]
        
        # Enflechamento da primeira seção
        if dat.or_sweep1 == 'Z':
            pop[i].sweep1 = 'Z'
        else:
            pop[i].sweep1 = dat.or_sweep1 + random.random()*dat.sweep1_ext[random.randint(0,1)]
        
        # Enflechamento da segunda seção
        if dat.or_sweep2 == 'Z':
            pop[i].sweep2 = 'Z'
        else:
            pop[i].sweep2 = dat.or_sweep2 + random.random()*dat.sweep2_ext[random.randint(0,1)]
    
    else: # Estabelecer enflechamento da asa trapezoidal simples 
        if dat.or_sweep == 'Z':
            pop[i].sweep = 'Z'
        else:
            pop[i].sweep = dat.or_sweep + random.random()*dat.sweep_ext[random.randint(0,1)]
        
# Gerar struct que guarda o melhor perfil de cada geração        
archive = empty.repeat(dat.iter)

## Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for loop in range(dat.iter):
    
    # Simular as asas e obter dados
    select = np.ones((1,dat.N))
    print('<< Simulação das asas >>')
    if dat.method == 1: # Painéis 3D
        for i in range(dat.N):
            print('Indivíduo ' + str(i+1))
            pop[i] = run_apame_mesher_cst(pop[i],dat,1) # Obter malha do apame
            pop[i].aero = run_apame(pop[i],dat) # Fazer simulação
            
            # Marcar indivíduos que não tenham convergido na simulação
            if isinstance(pop[i].aero,str) and pop[i].aero == 'n':
                select[[0],[i]] = 0
            
    elif dat.method == 2: # VLM
        for i in range(dat.N):
            print('Indivíduo ' + str(i+1))
            pop[i] = run_openvsp_cst_VLM_specify(pop[i],dat) # Obter geometria OpenVSP            
            pop[i].aero = run_vspaero_VLM(pop[i],dat) # Simulação
            
            # Marcar indivíduos que não tenham convergido na simulação
            if isinstance(pop[i].aero,str) and pop[i].aero == 'n':
                select[[0],[i]] = 0
    
    # Encontrar indivíduos problemáticos para já ignorá-los durante a atribuição
    # de pontuação
    # select contém os indivíduos com aero = 'n'; select2 contém os outros, que
    # convergiram nas simulações
    select = np.argwhere(select==0)[:,1]
    select2 = np.arange(0,dat.N)
    select2 = np.delete(select2,select)

    if len(select) == dat.N:
        raise TypeError('Nenhum aerofólio convergiu nas simulações')
    
     # Atribuir pontuações (fitnesses)
    pop = fitness_cst_3D(pop,dat,select2);
    
    # Isto serve pra pôr todas as pontuações em um vetor
    weights = np.zeros(dat.N)
    for i in range(dat.N):
        weights[i] = pop[i].score
    
    # Guardar a melhor asa de cada iteração
    pos = np.argmax(weights)
    archive[loop] = pop[pos].deepcopy()
    
    # Mostrar a melhor asa
    plt.figure()
    plot_planform(pop[pos],'k')
    plt.axis('equal');plt.grid('true')
    plt.title('Iteração ' + str(loop+1))   

    # Parar o código aqui na última iteração, já que nesse cenário o resto
    # do código é inútil
    if loop == dat.iter: break

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
        
        
        # if pop[par[0]].type != pop[par[1]].type: # Caso as asas sejam de tipos diferentes
        #     op = random.randint(1,2)
        # else: # Caso as asas sejam de tipos iguais
        #     op = 1        
        
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
            if cross_op_v[1] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1: #&& pop[par[1]].b1 <= pop[par[2]].b - dat.b1_ext[3]*2 && pop[par[2]].b1 <= pop[par[1]].b - dat.b1_ext[3]*2
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
            if cross_op_v[4] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1 and not dat.or_c_m == 'L':
                chi[c].c_m = pop[par[1]].c_m
                chi[c+1].c_m = pop[par[0]].c_m
            
            # Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi[c].type == 1:
                
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[c].b1 > chi[c].b-dat.b1_ext[3]*2:
                    chi[c].b1 = chi[c].b-dat.b1_ext[3]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[c].c_m > chi[c].c_r and dat.or_c_m != 'L':
                    chi[c].c_m = chi[c].c_r
            
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[c].c_m < chi[c].c_t and dat.or_c_m != 'L':
                    chi[c].c_m = chi[c].c_t
                
            if chi[c+1].type == 1:
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[c+1].b1 > chi[c+1].b-dat.b1_ext[3]*2:
                    chi[c+1].b1 = chi[c+1].b-dat.b1_ext[3]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[c+1].c_m > chi[c+1].c_r and dat.or_c_m != 'L':
                    chi[c+1].c_m = chi[c+1].c_r
                
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[c+1].c_m < chi[c+1].c_t and dat.or_c_m != 'L':
                    chi[c+1].c_m = chi[c+1].c_t
            
            # Aerofólio da raiz
            if cross_op_v[5] == 1:
                af1_ex = chi[c].v_ex_r
                af1_in = chi[c].v_in_r
                af2_ex = chi[c+1].v_ex_r
                af2_in = chi[c+1].v_in_r
                
                check = 0
                while check == 0:
                
                    m = random.randint(1,2)
                    
                    if m == 1: # Trocar os perfis completamente
                        chi[c].v_ex_r = af2_ex
                        chi[c].v_r_in = af2_in
                        chi[c+1].v_ex_r = af1_ex
                        chi[c+1].v_r_in = af1_in
                    
                    elif m == 2: # Fazer o crossover das características
                            
                        # v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                        if chi[c].symm_r == 0 and chi[c+1].symm_r == 0: # Se ambos forem assimétricos
                            n = random.randint(1,4)

                            if n == 1: # Trocar os extradorsos e intradorsos inteiros
                                temp1_ex = af1_ex
                                temp1_in = af2_in
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                                
                            elif n == 2: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = np.hstack((af1_in[0],af2_in[1:]))
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = np.hstack((af2_in[0],af1_in[1:]))
                                
                            elif n == 3: # Trocar os pesos intermediários
                                if dat.BPn_r == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_r],af2_ex[dat.BPn_r:]))
                                    temp1_in = np.hstack((af2_in[0],af1_in[1:dat.BPn_r],af2_in[dat.BPn_r:]))
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_r],af1_ex[dat.BPn_r:]))
                                    temp2_in = np.hstack((af1_in[0],af2_in[1:dat.BPn_r],af1_in[dat.BPn_r:]))
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(1,dat.BPn_r-1)
                                    num2 = random.randint(1,dat.BPn_r-1)
                                    temp1_1 = np.hstack((af1_ex[1:num1+1],af2_ex[num1+1:dat.BPn_r]))
                                    temp1_2 = np.hstack((af1_in[1:num2+1],af2_in[num2+1:dat.BPn_r]))
                                    temp2_1 = np.hstack((af2_ex[1:num1+1],af1_ex[num1+1:dat.BPn_r]))
                                    temp2_2 = np.hstack((af2_in[1:num2+1],af1_in[num2+1:dat.BPn_r]))
                                    temp1_ex = np.hstack(([af1_ex[0],temp2_1,af1_ex[dat.BPn_r:]]))
                                    temp1_in = np.hstack(([af1_in[0],temp2_2,af1_in[dat.BPn_r:]]))
                                    temp2_ex = np.hstack(([af2_ex[0],temp1_1,af2_ex[dat.BPn_r:]]))
                                    temp2_in = np.hstack(([af2_in[0],temp1_2,af2_in[dat.BPn_r:]]))
                                
                            elif n == 4: # Trocar os ângulos do bordo de fuga
                            
                                temp1_ex = np.hstack(([af2_ex[0:dat.BPn_r],af1_ex[dat.BPn_r],af2_ex[dat.BPn_r+1:]]))
                                temp1_in = np.hstack(([af2_in[0:dat.BPn_r],af1_in[dat.BPn_r],af2_in[dat.BPn_r+1:]]))
                                temp2_ex = np.hstack(([af1_ex[0:dat.BPn_r],af2_ex[dat.BPn_r],af1_ex[dat.BPn_r+1:]]))
                                temp2_in = np.hstack(([af1_in[0:dat.BPn_r],af2_in[dat.BPn_r],af1_in[dat.BPn_r+1:]]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_r = temp1_ex
                                chi[c].v_in_r = temp1_in
                                chi[c+1].v_ex_r = temp2_ex
                                chi[c+1].v_in_r = temp2_in
                            else:
                                chi[c].v_ex_r = temp2_ex
                                chi[c].v_in_r = temp2_in
                                chi[c+1].v_ex_r = temp1_ex
                                chi[c+1].v_in_r = temp1_ex
                            
                            if n == 1:
                                # Checar a separação dos bordos de fuga. Se não cumprirem o 
                                # requisito de separação, alterar o ângulo do bordo de fuga
                                # do intradorso
                                if chi[c].v_in_r[dat.BPn_r] < [dat.B_ext2_r[0]-chi[c].v_ex_r[dat.BPn_r]]:
                                    chi[c].v_in_r[dat.BPn_r] = dat.B_ext2_r[0]-chi[c].v_ex_r[dat.BPn_r]
                                
                                if chi[c+1].v_in_r[dat.BPn_r] < [dat.B_ext2_r[0]-chi[c+1].v_ex_r[dat.BPn_r]]:
                                    chi[c+1].v_in_r[dat.BPn_r] = dat.B_ext2_r[0]-chi[c+1].v_ex_r[dat.BPn_r]
                                
                                # Decisão de alterar o intradorso em base de uma nota na página
                                # 57(87) do Raymer (2018)
                            
                            
                            # Checar os pesos
                            sum1 = sum(chi[c].v_ex_r[1:dat.BPn_r])
                            sum2 = sum(chi[c].v_in_r[1:dat.BPn_r])
                            if sum2 > sum1:continue
                            sum1 = sum(chi[c+1].v_ex_r[1:dat.BPn_r])
                            sum2 = sum(chi[c+1].v_in_r[1:dat.BPn_r])
                            if sum2 > sum1:continue
                            
                            # Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                            # for usado como variável)
                            #chi(c).v(dat.BPn_r) = -chi(c).v_ex(dat.BPn_r)*chi(c).v_ex(dat.BPn_r-1)/chi(c).v(dat.BPn_r-1);
                            #chi(c+1).v(dat.BPn_r) = -chi(c+1).v_ex(dat.BPn_r)*chi(c+1).v_ex(dat.BPn_r-1)/chi(c+1).v(dat.BPn_r-1);
                        
                            chi[c].symm_r = 0
                            chi[c+1].symm_r = 0
                    
                        elif chi[c].symm_r == 1 and chi[c+1].symm_r == 1: # Se ambos forem simétricos
                            n = random.randint(1,3)
                                          
                            if n == 1: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = chi[c].v_ex_r
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = chi[c+1].v_ex_r    
                            
                            elif n == 2: # Trocar os pesos intermediários
                            
                                if dat.BPn_r == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_r],af2_ex[dat.BPn_r:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_r],af1_ex[dat.BPn_r:]))
                                    temp2_in = temp2_ex
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(2,dat.BPn_r) 
                                    temp1_1 = np.hstack((af1_ex[1:num1],af2_ex[num1:dat.BPn_r]))
                                    temp2_1 = np.hstack((af2_ex[1:num1],af1_ex[num1:dat.BPn_r]))
                                    temp1_ex = np.hstack((af1_ex[0],temp2_1,af1_ex[dat.BPn_r:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af2_ex[0],temp1_1,af2_ex[dat.BPn_r:]))
                                    temp2_in = temp2_ex
                                
                            if n == 3: # Trocar os ângulos do bordo de fuga
                                temp1_ex = np.hstack((af2_ex[0:dat.BPn_r],af1_ex[dat.BPn_r],af2_ex[dat.BPn_r+1]))
                                temp1_in = np.hstack((af2_in[0:dat.BPn_r],af1_in[dat.BPn_r],af2_in[dat.BPn_r+1]))
                                temp2_ex = np.hstack((af1_ex[0:dat.BPn_r],af2_ex[dat.BPn_r],af1_ex[dat.BPn_r+1]))
                                temp2_in = np.hstack((af1_in[0:dat.BPn_r],af2_in[dat.BPn_r],af1_in[dat.BPn_r+1]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_r = temp1_ex
                                chi[c].v_in_r = temp1_in
                                chi[c+1].v_ex_r = temp2_ex
                                chi[c+1].v_in_r = temp2_in
                            else:
                                chi[c].v_ex_r = temp2_ex
                                chi[c].v_in_r = temp2_in
                                chi[c+1].v_ex_r = temp1_ex
                                chi[c+1].v_in_r = temp1_ex
                            
                            chi[c].symm = 1
                            chi[c+1].symm = 1
                        
                        else: # Se um for simétrico e o outro for assimétrico
                            # Transformar o simétrico em um assimétrico e vice-versa
                            if pop[par[0]].symm == 0:
                                temp1_ex = af1_ex
                                temp1_in = af1_ex
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                            else:
                                temp1_ex = af2_ex
                                temp1_in = af2_ex
                                temp2_ex = af1_ex
                                temp2_in = af2_in

                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_r = temp1_ex
                                chi[c].v_in_r = temp1_in
                                chi[c+1].v_ex_r = temp2_ex
                                chi[c+1].v_in_r = temp2_in
                            else:
                                chi[c].v_ex_r = temp2_ex
                                chi[c].v_in_r = temp2_in
                                chi[c+1].v_ex_r = temp1_ex
                                chi[c+1].v_in_r = temp1_in
                            
                            # Checar a separação dos bordos de fuga. Se não cumprirem o 
                            # requisito de separação, alterar o ângulo do bordo de fuga
                            # do intradorso
                            if chi[c].v_in_r[dat.BPn_r] < [dat.B_ext2_r[0]-chi[c].v_ex_r[dat.BPn_r]]:
                                chi[c].v_in_r[dat.BPn_r] = dat.B_ext2_r[0]-chi[c].v_ex_r[dat.BPn_r]
                            
                            if chi[c+1].v_in_r[dat.BPn_r] < [dat.B_ext2_r[0]-chi[c+1].v_ex_r[dat.BPn_r]]:
                                chi[c+1].v_in_r[dat.BPn_r] = dat.B_ext2_r[0]-chi[c+1].v_ex_r[dat.BPn_r]
                            
                            # Decisão de alterar o intradorso em base de uma nota na página
                            # 57(87) do Raymer (2018)
                            
                            # Checar os pesos
                            sum1 = sum(chi[c+1].v_ex_r[1:dat.BPn_r])
                            sum2 = sum(chi[c+1].v_in_r[1:dat.BPn_r])
                            if sum2 > sum1:continue
                            
                            chi[c].symm_r = 1
                            chi[c+1].symm_r = 0 
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi[c].v_ex_r,chi[c].v_in_r,dat,[dat.N1_r,dat.N2_r]),dat)
                    if check == 0:continue
                    check = quality(run_cst_TCC2_3D(chi[c+1].v_ex_r,chi[c+1].v_in_r,dat,[dat.N1_r,dat.N2_r]),dat)
                
           
            
            # Aerofólio do meio (asas bitrapezoidais apenas) (ignorar caso a opção 'L' seja aplicada ao aerofólio do meio)
            if cross_op_v[6] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1 and not 'L' in dat.le_R_ext1_m:
                af1_ex = chi[c].v_ex_m
                af1_in = chi[c].v_in_m
                af2_ex = chi[c+1].v_ex_m
                af2_in = chi[c+1].v_in_m
                
                check = 0
                while check == 0:
                
                    m = random.randint(1,2)
                    
                    if m == 1: # Trocar os perfis completamente
                        chi[c].v_ex_m = af2_ex
                        chi[c].v_m_in = af2_in
                        chi[c+1].v_ex_m = af1_ex
                        chi[c+1].v_m_in = af1_in
                    
                    elif m == 2: # Fazer o crossover das características
                            
                        # v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                        if chi[c].symm_m == 0 and chi[c+1].symm_m == 0: # Se ambos forem assimétricos
                            n = random.randint(1,4)

                            if n == 1: # Trocar os extradorsos e intradorsos inteiros
                                temp1_ex = af1_ex
                                temp1_in = af2_in
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                                
                            elif n == 2: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = np.hstack((af1_in[0],af2_in[1:]))
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = np.hstack((af2_in[0],af1_in[1:]))
                                
                            elif n == 3: # Trocar os pesos intermediários
                                if dat.BPn_m == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_m],af2_ex[dat.BPn_m:]))
                                    temp1_in = np.hstack((af2_in[0],af1_in[1:dat.BPn_m],af2_in[dat.BPn_m:]))
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_m],af1_ex[dat.BPn_m:]))
                                    temp2_in = np.hstack((af1_in[0],af2_in[1:dat.BPn_m],af1_in[dat.BPn_m:]))
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(1,dat.BPn_m-1)
                                    num2 = random.randint(1,dat.BPn_m-1)
                                    temp1_1 = np.hstack((af1_ex[1:num1+1],af2_ex[num1+1:dat.BPn_m]))
                                    temp1_2 = np.hstack((af1_in[1:num2+1],af2_in[num2+1:dat.BPn_m]))
                                    temp2_1 = np.hstack((af2_ex[1:num1+1],af1_ex[num1+1:dat.BPn_m]))
                                    temp2_2 = np.hstack((af2_in[1:num2+1],af1_in[num2+1:dat.BPn_m]))
                                    temp1_ex = np.hstack(([af1_ex[0],temp2_1,af1_ex[dat.BPn_m:]]))
                                    temp1_in = np.hstack(([af1_in[0],temp2_2,af1_in[dat.BPn_m:]]))
                                    temp2_ex = np.hstack(([af2_ex[0],temp1_1,af2_ex[dat.BPn_m:]]))
                                    temp2_in = np.hstack(([af2_in[0],temp1_2,af2_in[dat.BPn_m:]]))
                                
                            elif n == 4: # Trocar os ângulos do bordo de fuga
                            
                                temp1_ex = np.hstack(([af2_ex[0:dat.BPn_m],af1_ex[dat.BPn_m],af2_ex[dat.BPn_m+1:]]))
                                temp1_in = np.hstack(([af2_in[0:dat.BPn_m],af1_in[dat.BPn_m],af2_in[dat.BPn_m+1:]]))
                                temp2_ex = np.hstack(([af1_ex[0:dat.BPn_m],af2_ex[dat.BPn_m],af1_ex[dat.BPn_m+1:]]))
                                temp2_in = np.hstack(([af1_in[0:dat.BPn_m],af2_in[dat.BPn_m],af1_in[dat.BPn_m+1:]]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_m = temp1_ex
                                chi[c].v_in_m = temp1_in
                                chi[c+1].v_ex_m = temp2_ex
                                chi[c+1].v_in_m = temp2_in
                            else:
                                chi[c].v_ex_m = temp2_ex
                                chi[c].v_in_m = temp2_in
                                chi[c+1].v_ex_m = temp1_ex
                                chi[c+1].v_in_m = temp1_ex
                            
                            if n == 1:
                                # Checar a separação dos bordos de fuga. Se não cumprirem o 
                                # requisito de separação, alterar o ângulo do bordo de fuga
                                # do intradorso
                                if chi[c].v_in_m[dat.BPn_m] < [dat.B_ext2_m[0]-chi[c].v_ex_m[dat.BPn_m]]:
                                    chi[c].v_in_m[dat.BPn_m] = dat.B_ext2_m[0]-chi[c].v_ex_m[dat.BPn_m]
                                
                                if chi[c+1].v_in_m[dat.BPn_m] < [dat.B_ext2_m[0]-chi[c+1].v_ex_m[dat.BPn_m]]:
                                    chi[c+1].v_in_m[dat.BPn_m] = dat.B_ext2_m[0]-chi[c+1].v_ex_m[dat.BPn_m]
                                
                                # Decisão de alterar o intradorso em base de uma nota na página
                                # 57(87) do Raymer (2018)
                            
                            
                            # Checar os pesos
                            sum1 = sum(chi[c].v_ex_m[1:dat.BPn_m])
                            sum2 = sum(chi[c].v_in_m[1:dat.BPn_m])
                            if sum2 > sum1:continue
                            sum1 = sum(chi[c+1].v_ex_m[1:dat.BPn_m])
                            sum2 = sum(chi[c+1].v_in_m[1:dat.BPn_m])
                            if sum2 > sum1:continue
                            
                            # Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                            # for usado como variável)
                            #chi(c).v(dat.BPn_m) = -chi(c).v_ex(dat.BPn_m)*chi(c).v_ex(dat.BPn_m-1)/chi(c).v(dat.BPn_m-1);
                            #chi(c+1).v(dat.BPn_m) = -chi(c+1).v_ex(dat.BPn_m)*chi(c+1).v_ex(dat.BPn_m-1)/chi(c+1).v(dat.BPn_m-1);
                        
                            chi[c].symm_m = 0
                            chi[c+1].symm_m = 0
                    
                        elif chi[c].symm_m == 1 and chi[c+1].symm_m == 1: # Se ambos forem simétricos
                            n = random.randint(1,3)
                                          
                            if n == 1: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = chi[c].v_ex_m
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = chi[c+1].v_ex_m    
                            
                            elif n == 2: # Trocar os pesos intermediários
                            
                                if dat.BPn_m == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_m],af2_ex[dat.BPn_m:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_m],af1_ex[dat.BPn_m:]))
                                    temp2_in = temp2_ex
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(2,dat.BPn_m) 
                                    temp1_1 = np.hstack((af1_ex[1:num1],af2_ex[num1:dat.BPn_m]))
                                    temp2_1 = np.hstack((af2_ex[1:num1],af1_ex[num1:dat.BPn_m]))
                                    temp1_ex = np.hstack((af1_ex[0],temp2_1,af1_ex[dat.BPn_m:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af2_ex[0],temp1_1,af2_ex[dat.BPn_m:]))
                                    temp2_in = temp2_ex
                                
                            if n == 3: # Trocar os ângulos do bordo de fuga
                                temp1_ex = np.hstack((af2_ex[0:dat.BPn_m],af1_ex[dat.BPn_m],af2_ex[dat.BPn_m+1]))
                                temp1_in = np.hstack((af2_in[0:dat.BPn_m],af1_in[dat.BPn_m],af2_in[dat.BPn_m+1]))
                                temp2_ex = np.hstack((af1_ex[0:dat.BPn_m],af2_ex[dat.BPn_m],af1_ex[dat.BPn_m+1]))
                                temp2_in = np.hstack((af1_in[0:dat.BPn_m],af2_in[dat.BPn_m],af1_in[dat.BPn_m+1]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_m = temp1_ex
                                chi[c].v_in_m = temp1_in
                                chi[c+1].v_ex_m = temp2_ex
                                chi[c+1].v_in_m = temp2_in
                            else:
                                chi[c].v_ex_m = temp2_ex
                                chi[c].v_in_m = temp2_in
                                chi[c+1].v_ex_m = temp1_ex
                                chi[c+1].v_in_m = temp1_ex
                            
                            chi[c].symm = 1
                            chi[c+1].symm = 1
                        
                        else: # Se um for simétrico e o outro for assimétrico
                            # Transformar o simétrico em um assimétrico e vice-versa
                            if pop[par[0]].symm == 0:
                                temp1_ex = af1_ex
                                temp1_in = af1_ex
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                            else:
                                temp1_ex = af2_ex
                                temp1_in = af2_ex
                                temp2_ex = af1_ex
                                temp2_in = af2_in

                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_m = temp1_ex
                                chi[c].v_in_m = temp1_in
                                chi[c+1].v_ex_m = temp2_ex
                                chi[c+1].v_in_m = temp2_in
                            else:
                                chi[c].v_ex_m = temp2_ex
                                chi[c].v_in_m = temp2_in
                                chi[c+1].v_ex_m = temp1_ex
                                chi[c+1].v_in_m = temp1_in
                            
                            # Checar a separação dos bordos de fuga. Se não cumprirem o 
                            # requisito de separação, alterar o ângulo do bordo de fuga
                            # do intradorso
                            if chi[c].v_in_m[dat.BPn_m] < [dat.B_ext2_m[0]-chi[c].v_ex_m[dat.BPn_m]]:
                                chi[c].v_in_m[dat.BPn_m] = dat.B_ext2_m[0]-chi[c].v_ex_m[dat.BPn_m]
                            
                            if chi[c+1].v_in_m[dat.BPn_m] < [dat.B_ext2_m[0]-chi[c+1].v_ex_m[dat.BPn_m]]:
                                chi[c+1].v_in_m[dat.BPn_m] = dat.B_ext2_m[0]-chi[c+1].v_ex_m[dat.BPn_m]
                            
                            # Decisão de alterar o intradorso em base de uma nota na página
                            # 57(87) do Raymer (2018)
                            
                            # Checar os pesos
                            sum1 = sum(chi[c+1].v_ex_m[1:dat.BPn_m])
                            sum2 = sum(chi[c+1].v_in_m[1:dat.BPn_m])
                            if sum2 > sum1:continue
                            
                            chi[c].symm_m = 1
                            chi[c+1].symm_m = 0 
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi[c].v_ex_m,chi[c].v_in_m,dat,[dat.N1_m,dat.N2_m]),dat)
                    if check == 0:continue
                    check = quality(run_cst_TCC2_3D(chi[c+1].v_ex_m,chi[c+1].v_in_m,dat,[dat.N1_m,dat.N2_m]),dat)
            
            # Aerofólio da ponta
            if cross_op_v[7] == 1:
                af1_ex = chi[c].v_ex_t
                af1_in = chi[c].v_in_t
                af2_ex = chi[c+1].v_ex_t
                af2_in = chi[c+1].v_in_t
                
                check = 0
                while check == 0:
                
                    m = random.randint(1,2)
                    
                    if m == 1: # Trocar os perfis completamente
                        chi[c].v_ex_t = af2_ex
                        chi[c].v_t_in = af2_in
                        chi[c+1].v_ex_t = af1_ex
                        chi[c+1].v_t_in = af1_in
                    
                    elif m == 2: # Fazer o crossover das características
                            
                        # v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                        if chi[c].symm_t == 0 and chi[c+1].symm_t == 0: # Se ambos forem assimétricos
                            n = random.randint(1,4)

                            if n == 1: # Trocar os extradorsos e intradorsos inteiros
                                temp1_ex = af1_ex
                                temp1_in = af2_in
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                                
                            elif n == 2: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = np.hstack((af1_in[0],af2_in[1:]))
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = np.hstack((af2_in[0],af1_in[1:]))
                                
                            elif n == 3: # Trocar os pesos intermediários
                                if dat.BPn_t == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_t],af2_ex[dat.BPn_t:]))
                                    temp1_in = np.hstack((af2_in[0],af1_in[1:dat.BPn_t],af2_in[dat.BPn_t:]))
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_t],af1_ex[dat.BPn_t:]))
                                    temp2_in = np.hstack((af1_in[0],af2_in[1:dat.BPn_t],af1_in[dat.BPn_t:]))
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(1,dat.BPn_t-1)
                                    num2 = random.randint(1,dat.BPn_t-1)
                                    temp1_1 = np.hstack((af1_ex[1:num1+1],af2_ex[num1+1:dat.BPn_t]))
                                    temp1_2 = np.hstack((af1_in[1:num2+1],af2_in[num2+1:dat.BPn_t]))
                                    temp2_1 = np.hstack((af2_ex[1:num1+1],af1_ex[num1+1:dat.BPn_t]))
                                    temp2_2 = np.hstack((af2_in[1:num2+1],af1_in[num2+1:dat.BPn_t]))
                                    temp1_ex = np.hstack(([af1_ex[0],temp2_1,af1_ex[dat.BPn_t:]]))
                                    temp1_in = np.hstack(([af1_in[0],temp2_2,af1_in[dat.BPn_t:]]))
                                    temp2_ex = np.hstack(([af2_ex[0],temp1_1,af2_ex[dat.BPn_t:]]))
                                    temp2_in = np.hstack(([af2_in[0],temp1_2,af2_in[dat.BPn_t:]]))
                                
                            elif n == 4: # Trocar os ângulos do bordo de fuga
                            
                                temp1_ex = np.hstack(([af2_ex[0:dat.BPn_t],af1_ex[dat.BPn_t],af2_ex[dat.BPn_t+1:]]))
                                temp1_in = np.hstack(([af2_in[0:dat.BPn_t],af1_in[dat.BPn_t],af2_in[dat.BPn_t+1:]]))
                                temp2_ex = np.hstack(([af1_ex[0:dat.BPn_t],af2_ex[dat.BPn_t],af1_ex[dat.BPn_t+1:]]))
                                temp2_in = np.hstack(([af1_in[0:dat.BPn_t],af2_in[dat.BPn_t],af1_in[dat.BPn_t+1:]]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_t = temp1_ex
                                chi[c].v_in_t = temp1_in
                                chi[c+1].v_ex_t = temp2_ex
                                chi[c+1].v_in_t = temp2_in
                            else:
                                chi[c].v_ex_t = temp2_ex
                                chi[c].v_in_t = temp2_in
                                chi[c+1].v_ex_t = temp1_ex
                                chi[c+1].v_in_t = temp1_ex
                            
                            if n == 1:
                                # Checar a separação dos bordos de fuga. Se não cumprirem o 
                                # requisito de separação, alterar o ângulo do bordo de fuga
                                # do intradorso
                                if chi[c].v_in_t[dat.BPn_t] < [dat.B_ext2_t[0]-chi[c].v_ex_t[dat.BPn_t]]:
                                    chi[c].v_in_t[dat.BPn_t] = dat.B_ext2_t[0]-chi[c].v_ex_t[dat.BPn_t]
                                
                                if chi[c+1].v_in_t[dat.BPn_t] < [dat.B_ext2_t[0]-chi[c+1].v_ex_t[dat.BPn_t]]:
                                    chi[c+1].v_in_t[dat.BPn_t] = dat.B_ext2_t[0]-chi[c+1].v_ex_t[dat.BPn_t]
                                
                                # Decisão de alterar o intradorso em base de uma nota na página
                                # 57(87) do Raymer (2018)
                            
                            
                            # Checar os pesos
                            sum1 = sum(chi[c].v_ex_t[1:dat.BPn_t])
                            sum2 = sum(chi[c].v_in_t[1:dat.BPn_t])
                            if sum2 > sum1:continue
                            sum1 = sum(chi[c+1].v_ex_t[1:dat.BPn_t])
                            sum2 = sum(chi[c+1].v_in_t[1:dat.BPn_t])
                            if sum2 > sum1:continue
                            
                            # Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                            # for usado como variável)
                            #chi(c).v(dat.BPn_t) = -chi(c).v_ex(dat.BPn_t)*chi(c).v_ex(dat.BPn_t-1)/chi(c).v(dat.BPn_t-1);
                            #chi(c+1).v(dat.BPn_t) = -chi(c+1).v_ex(dat.BPn_t)*chi(c+1).v_ex(dat.BPn_t-1)/chi(c+1).v(dat.BPn_t-1);
                        
                            chi[c].symm_t = 0
                            chi[c+1].symm_t = 0
                    
                        elif chi[c].symm_t == 1 and chi[c+1].symm_t == 1: # Se ambos forem simétricos
                            n = random.randint(1,3)
                                          
                            if n == 1: # Trocar o raio do bordo de ataque
                                temp1_ex = np.hstack((af1_ex[0],af2_ex[1:]))
                                temp1_in = chi[c].v_ex_t
                                temp2_ex = np.hstack((af2_ex[0],af1_ex[1:]))
                                temp2_in = chi[c+1].v_ex_t    
                            
                            elif n == 2: # Trocar os pesos intermediários
                            
                                if dat.BPn_t == 2:
                                    op = 1
                                else:
                                    op = random.randint(1,2)
                                
                                if op == 1: # Trocar tudo    
                                    temp1_ex = np.hstack((af2_ex[0],af1_ex[1:dat.BPn_t],af2_ex[dat.BPn_t:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af1_ex[0],af2_ex[1:dat.BPn_t],af1_ex[dat.BPn_t:]))
                                    temp2_in = temp2_ex
                                    
                                else: # Trocar cortes
                                    num1 = random.randint(2,dat.BPn_t) 
                                    temp1_1 = np.hstack((af1_ex[1:num1],af2_ex[num1:dat.BPn_t]))
                                    temp2_1 = np.hstack((af2_ex[1:num1],af1_ex[num1:dat.BPn_t]))
                                    temp1_ex = np.hstack((af1_ex[0],temp2_1,af1_ex[dat.BPn_t:]))
                                    temp1_in = temp1_ex
                                    temp2_ex = np.hstack((af2_ex[0],temp1_1,af2_ex[dat.BPn_t:]))
                                    temp2_in = temp2_ex
                                
                            if n == 3: # Trocar os ângulos do bordo de fuga
                                temp1_ex = np.hstack((af2_ex[0:dat.BPn_t],af1_ex[dat.BPn_t],af2_ex[dat.BPn_t+1]))
                                temp1_in = np.hstack((af2_in[0:dat.BPn_t],af1_in[dat.BPn_t],af2_in[dat.BPn_t+1]))
                                temp2_ex = np.hstack((af1_ex[0:dat.BPn_t],af2_ex[dat.BPn_t],af1_ex[dat.BPn_t+1]))
                                temp2_in = np.hstack((af1_in[0:dat.BPn_t],af2_in[dat.BPn_t],af1_in[dat.BPn_t+1]))
                            
                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_t = temp1_ex
                                chi[c].v_in_t = temp1_in
                                chi[c+1].v_ex_t = temp2_ex
                                chi[c+1].v_in_t = temp2_in
                            else:
                                chi[c].v_ex_t = temp2_ex
                                chi[c].v_in_t = temp2_in
                                chi[c+1].v_ex_t = temp1_ex
                                chi[c+1].v_in_t = temp1_ex
                            
                            chi[c].symm = 1
                            chi[c+1].symm = 1
                        
                        else: # Se um for simétrico e o outro for assimétrico
                            # Transformar o simétrico em um assimétrico e vice-versa
                            if pop[par[0]].symm == 0:
                                temp1_ex = af1_ex
                                temp1_in = af1_ex
                                temp2_ex = af2_ex
                                temp2_in = af1_in
                            else:
                                temp1_ex = af2_ex
                                temp1_in = af2_ex
                                temp2_ex = af1_ex
                                temp2_in = af2_in

                            # Decidir pra qual asa vai cada um dos aerofólios novos
                            if random.randint(0,1) == 1:
                                chi[c].v_ex_t = temp1_ex
                                chi[c].v_in_t = temp1_in
                                chi[c+1].v_ex_t = temp2_ex
                                chi[c+1].v_in_t = temp2_in
                            else:
                                chi[c].v_ex_t = temp2_ex
                                chi[c].v_in_t = temp2_in
                                chi[c+1].v_ex_t = temp1_ex
                                chi[c+1].v_in_t = temp1_in
                            
                            # Checar a separação dos bordos de fuga. Se não cumprirem o 
                            # requisito de separação, alterar o ângulo do bordo de fuga
                            # do intradorso
                            if chi[c].v_in_t[dat.BPn_t] < [dat.B_ext2_t[0]-chi[c].v_ex_t[dat.BPn_t]]:
                                chi[c].v_in_t[dat.BPn_t] = dat.B_ext2_t[0]-chi[c].v_ex_t[dat.BPn_t]
                            
                            if chi[c+1].v_in_t[dat.BPn_t] < [dat.B_ext2_t[0]-chi[c+1].v_ex_t[dat.BPn_t]]:
                                chi[c+1].v_in_t[dat.BPn_t] = dat.B_ext2_t[0]-chi[c+1].v_ex_t[dat.BPn_t]
                            
                            # Decisão de alterar o intradorso em base de uma nota na página
                            # 57(87) do Raymer (2018)
                            
                            # Checar os pesos
                            sum1 = sum(chi[c+1].v_ex_t[1:dat.BPn_t])
                            sum2 = sum(chi[c+1].v_in_t[1:dat.BPn_t])
                            if sum2 > sum1:continue
                            
                            chi[c].symm_t = 1
                            chi[c+1].symm_t = 0 
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi[c].v_ex_t,chi[c].v_in_t,dat,[dat.N1_t,dat.N2_t]),dat)
                    if check == 0:continue
                    check = quality(run_cst_TCC2_3D(chi[c+1].v_ex_t,chi[c+1].v_in_t,dat,[dat.N1_t,dat.N2_t]),dat)     
            
            # Torção geométrica na ponta
            if cross_op_v[8] == 1: 
                chi[c].tw_t = pop[par[1]].tw_t
                chi[c+1].tw_t = pop[par[0]].tw_t       
                
            # Enflechamento total (asas trapezoidais simples)
            if cross_op_v[9] == 1 and pop[par[0]].type == 0 and pop[par[1]].type == 0:
                chi[c].sweep = pop[par[1]].sweep;
                chi[c+1].sweep = pop[par[0]].sweep;
        
            # Enflechamento da primeira seção (asas bitrapezoidais)
            if cross_op_v[10] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1:
                chi[c].sweep1 = pop[par[1]].sweep1;
                chi[c+1].sweep1 = pop[par[0]].sweep1;
            
            # Enflechamento da segunda seção (asas bitrapezoidais)
            if cross_op_v[11] == 1 and pop[par[0]].type == 1 and pop[par[1]].type == 1:
                chi[c].sweep2 = pop[par[1]].sweep2;
                chi[c+1].sweep2 = pop[par[0]].sweep2;
            
        # else: # Transformar uma asa trapezoidal simples em trapezoidal dupla e vice-versa
        #     # (adiciona-se/retira-se o aerofólio do meio)
        #     if pop[par[0]].type == 0: # Se a primeira do par for trapezoidal simples
        #         chi[c] = pop[par[0]]
        #         chi[c].type = 1
        #         chi[c].b1 = pop[par[1]].b1
        #         chi[c].c_m = pop[par[1]].c_m
        #         chi[c].af_m = pop[par[1]].af_m
        #         chi[c].v_ex_m = pop[par[1]].v_ex_m
        #         chi[c].v_m = pop[par[1]].v_m
        #         chi[c].tw_m = pop[par[1]].tw_m
        #         chi[c].sweep1 = pop[par[1]].sweep1
        #         chi[c].sweep2 = pop[par[1]].sweep2
                
        #         chi[c+1] = pop[par[1]]
        #         chi[c+1].type = 0
        #         chi[c+1].b1 = None
        #         chi[c+1].c_m = None
        #         chi[c+1].v_ex_m = np.zeros(1,dat.BPn_m+2)
        #         chi[c+1].v_m = np.zeros(1,dat.BPn_m+2)
        #         chi[c+1].tw_m = None
        #         chi[c+1].sweep = pop[par[1]].sweep
                
        #     else: # Se a primeira do par for trapezoidal dupla
        #         chi[c] = pop[par[1]]
        #         chi[c].type = 1
        #         chi[c].b1 = pop[par[0]].b1
        #         chi[c].c_m = pop[par[0]].c_m
        #         chi[c].af_m = pop[par[0]].af_m
        #         chi[c].tw_m = pop[par[0]].tw_m
        #         chi[c].sweep1 = pop[par[0]].sweep1
        #         chi[c].sweep1 = pop[par[0]].sweep1
                
        #         chi[c+1] = pop[par[0]]
        #         chi[c+1].type = 0
        #         chi[c+1].b1 = None
        #         chi[c+1].c_m = None
        #         chi[c+1].af_m = np.zeros(3)
        #         chi[c+1].tw_m = None
        #         chi[c+1].sweep = pop[par[1]].sweep
        
        
            
            
        #     # Checar o requisito das envergaduras
        #     chi[c].b1 = dat.or_b1 + random.random()*dat.b1_ext(random.randint(0,1))
        #     if chi[c].b1 < dat.b1_ext[2]:
        #         chi[c].b1 = dat.b1_ext[2]
            
        #     if chi[c].b1 > chi[c].b - 2*dat.b1_ext[3]:
        #         chi[c].b1 = chi[c].b - 2*dat.b1_ext[3]
            
        #     # Checar o requisito dos comprimentos de corda
        #     if not dat.or_c_m == 'L':
        #         if chi[c].c_t > chi[c].c_m or chi[c].c_m > chi[c].c_r:
        #             chi[c].c_m = dat.or_c_m + random.random()*dat.c_m_ext[random.randint(0,1)]
        #             if chi[c].c_m < dat.c_m_ext[2]:
        #                 chi[c].c_m = dat.c_m_ext[2]
                    
        #             if chi[c].c_m > chi[c].c_r:
        #                 chi[c].c_m = chi[c].c_r
                    
        #             if chi[c].c_t > chi[c].c_m:
        #                 chi[c].c_t = chi[c].c_m
                    
        
        # Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
        if chi[c].type == 1:
            if chi[c].b1 > chi[c].b-2*dat.b1_ext[3] or chi[c].c_t > chi[c].c_m or chi[c].c_m > chi[c].c_r:
                raise TypeError('Problema na geometria da planta [c]')
            
        if chi[c+1].type == 1:
            if chi[c+1].b1 > chi[c+1].b-2*dat.b1_ext[3] or chi[c+1].c_t > chi[c+1].c_m or chi[c+1].c_m > chi[c+1].c_r:
                raise TypeError('Problema na geometria da planta [c+1]')
            
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
                chi[s].b = dat.or_b + random.randint()*dat.b_ext(random.randint(0,1))
                if chi[s].b < dat.b_ext[2]:
                    chi[s].b = dat.b_ext[2]
            
            # Envergadura da primeira seção b1 (asas bitrapezoidais apenas)
            if mu_op_v[1] == 1 and chi[s].type == 1:
                chi[s].b1 = dat.or_b1 + random.random()*dat.b1_ext[random.randint(0,1)]
                if chi[s].b1 < dat.b1_ext[2]:
                    chi[s].b1 = dat.b1_ext[2]
            
            # Corda da raiz c_r
            if mu_op_v[2] == 1:
                chi[s].c_r = dat.or_c_r + random.random()*dat.c_r_ext(random.randint(0,1))
                if chi[s].c_r < dat.c_r_ext[2]:
                    chi[s].c_r = dat.c_r_ext[2]
            
            # Corda da ponta c_t (deve cumprir o requisito c_t <= c_r)
            if mu_op_v[3] == 1:
                chi[s].c_t = dat.or_c_t + random.random()*dat.c_t_ext[random.randint(0,1)]
                if chi[s].c_t < dat.c_t_ext[2]:
                    chi[s].c_t = dat.c_t_ext[2]
            
            # Corda do meio c_m (asas bitrapezoidais apenas)
            if mu_op_v[4] == 1 and not dat.or_c_m == 'L' and dat.or_c_m != 'L':
                chi[s].c_m = dat.or_c_m + random.random()*dat.c_m_ext[random.randint(0,1)]
                if chi[s].c_m < dat.c_m_ext[2]:
                    chi[s].c_m = dat.c_m_ext[2]
            
            # Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi[s].type == 1:
                
                # Se a envergadura da primeira seção for maior do que permitido
                # pelo requisito de separação, atribuir o máximo valor que 
                # cumpre o requisito
                if chi[s].b1 > chi[s].b-dat.b1_ext_in[2]*2:
                    chi[s].b1 = chi[s].b-dat.b1_ext_in[2]*2
                
                # Se a corda do meio for maior do que a corda da raiz, atribuir
                # o valor da raiz ao meio
                if chi[s].c_m > chi[s].c_r:
                    chi[s].c_m = chi[s].c_r
                
                # Se a corda do meio for menor que a corda da ponta, atribuir
                # o valor da ponta ao meio
                if chi[s].c_m < chi[s].c_t:
                    chi[s].c_m = chi[s].c_t
            
            # Aerofólio da raiz
            if mu_op_v[5] == 1:
                check = 0
                while check == 0:
                    
                    temp = chi[s].deepcopy()
                    if temp.symm_r == 1: # Caso o perfil seja simétrico
                        n =([1,4,5])[random.randint(0,2)]
                    else: # Caso o perfil seja assimétrico
                        n = random.randint(1,5)
                    
                    if n == 1: # Alterar o raio do bordo de ataque
                        if temp.symm_r == 1:
                            P = 1
                        else:
                            P = random.randint(1,4)
                        
                        if P == 1: # Mudar ambos para o mesmo valor
                            temp.v_ex_r[0] = dat.or_v_ex_r[0] + random.random()*dat.le_R_ext1_r[random.randint(0,1)]
                            temp.v_in_r[0] = temp.v_ex_r[0]
                        elif P == 2: # Mudar ambos independentemente 
                            temp.v_ex_r[0] = dat.or_v_ex_r[0] + random.random()*dat.le_R_ext1_r[random.randint(0,1)]
                            temp.v_in_r[0] = dat.or_v_in_r[0] + random.random()*dat.le_R_ext2_r[random.randint(0,1)]
                        elif P == 3: # Mudar do extradorso
                            temp.v_ex_r[0] = dat.or_v_ex_r[0] + random.random()*dat.le_R_ext1_r[random.randint(0,1)]
                        else: # Mudar do intradorso
                            temp.v_in_r[0] = dat.or_v_in_r[0] + random.random()*dat.le_R_ext2_r[random.randint(0,1)]
                            
                        # Checar o raio do bordo de ataque (simétricos)
                        if temp.symm_r == 1 and temp.v_ex_r[0] < dat.le_R_ext1_r[2]:
                            temp.v_ex_r[0] = dat.le_R_ext1_r[2]
                            temp.v_in_r[0] = temp.v_ex_r[0]
                        
                        # Checar o raio do bordo de ataque (assimétricos)
                        if temp.symm_r == 0 and temp.v_ex_r[0] < dat.le_R_ext1_r[2]:
                            temp.v_ex_r[0] = dat.le_R_ext1_r[2]
                        if temp.symm_r == 0 and temp.v_in_r[0] < dat.le_R_ext2_r[2]:
                            temp.v_in_r[0] = dat.le_R_ext2_r[2]
                        
                    elif n == 2: # Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_r):
                            temp.v_ex_r[a] = dat.or_v_ex_r[a] + random.random()*dat.A_ext1_r[random.randint(0,1)]
                    
                    elif n == 3: # Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_r):
                            temp.v_in_r[a] = dat.or_v_in[a] + random.random()*dat.A_ext2_r[random.randint(0,1)]
    
                    elif n == 4: # Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                        if temp.symm_r == 1: # Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                            for a in range(1,dat.BPn_r):
                                temp.v_ex_r[a] = dat.or_v_ex_r[a] + random.random()*dat.A_ext1_r[random.randint(0,1)]
                                temp.v_in_r[a] = temp.v_ex_r[a]
                            
                        else: # Perfis assimétricos ficam com pesos intermediários distintos
                            for a in range (1,dat.BPn_r):
                                temp.v_ex_r[a] = dat.or_v_ex_r[a] + random.random()*dat.A_ext1_r[random.randint(0,1)]
                            
                            for a in range (1,dat.BPn_r):
                                temp.v_in_r[a] =dat.or_v_in[a] + random.random()*dat.A_ext2_r[random.randint(0,1)]
    
                    elif n == 5: # Alterar o ângulo do bordo de fuga
                        if temp.symm_r == 1: # Se for simétrico
                           temp.v_ex_r[dat.BPn_r] = dat.or_v_ex_r[dat.BPn_r] + random.random()*dat.B_ext1_r[random.randint(0,1)]
                           temp.v_in_r[dat.BPn_r] = temp.v_ex_r[dat.BPn_r]
                               
                           # Checar separação do bordo de fuga
                           if temp.v_ex_r[dat.BPn_r] + temp.v_in_r[dat.BPn_r] < dat.B_ext2_r[0]:
                               temp.v_ex_r[dat.BPn_r] = dat.B_ext2_r[0]/2
                               temp.v_in_r[dat.BPn_r] = dat.B_ext2_r[0]/2
                           
                        else: # Se for assimétrico
                            temp.v_ex_r[dat.BPn_r] = dat.or_v_ex_r[dat.BPn_r] + random.random()*dat.B_ext1_r[random.randint(0,1)]
                            temp.v_in_r[dat.BPn_r] = (dat.B_ext2_r[0] - temp.v_ex_r[dat.BPn_r]) + random.random()*dat.B_ext2_r[1]
                            
                            # Checar separação do bordo de fuga [beta2>=L-beta1]
                            if temp.v_in_r[dat.BPn_r] < dat.B_ext2_r[0] - temp.v_ex_r[dat.BPn_r]:
                                temp.v_in_r[dat.BPn_r] = dat.B_ext2_r[0] - temp.v_ex_r[dat.BPn_r]


                    # Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    # igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_r[1:dat.BPn_r])
                    sum2 = sum(temp.v_in_r[1:dat.BPn_r])
                    if sum2 > sum1:continue
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_r,temp.v_in_r,dat,[dat.N1_r,dat.N2_r]),dat)
            
            
                chi[s] = temp.deepcopy()
                
            
            
            # Aerofólio do meio (asas bitrapezoidais apenas)
            if mu_op_v[6] == 1 and chi[s].type == 1 and not 'L' in dat.le_r_ext1_m:
                check = 0
                while check == 0:
                    
                    temp = chi[s].deepcopy()
                    if temp.symm_m == 1: # Caso o perfil seja simétrico
                        n =([1,4,5])[random.randint(0,2)]
                    else: # Caso o perfil seja assimétrico
                        n = random.randint(1,5)
                    
                    if n == 1: # Alterar o raio do bordo de ataque
                        if temp.symm_m == 1:
                            P = 1
                        else:
                            P = random.randint(1,4)
                        
                        if P == 1: # Mudar ambos para o mesmo valor
                            temp.v_ex_m[0] = dat.or_v_ex_m[0] + random.random()*dat.le_R_ext1_m[random.randint(0,1)]
                            temp.v_in_m[0] = temp.v_ex_m[0]
                        elif P == 2: # Mudar ambos independentemente 
                            temp.v_ex_m[0] = dat.or_v_ex_m[0] + random.random()*dat.le_R_ext1_m[random.randint(0,1)]
                            temp.v_in_m[0] = dat.or_v_in_m[0] + random.random()*dat.le_R_ext2_m[random.randint(0,1)]
                        elif P == 3: # Mudar do extradorso
                            temp.v_ex_m[0] = dat.or_v_ex_m[0] + random.random()*dat.le_R_ext1_m[random.randint(0,1)]
                        else: # Mudar do intradorso
                            temp.v_in_m[0] = dat.or_v_in_m[0] + random.random()*dat.le_R_ext2_m[random.randint(0,1)]
                            
                        # Checar o raio do bordo de ataque (simétricos)
                        if temp.symm_m == 1 and temp.v_ex_m[0] < dat.le_R_ext1_m[2]:
                            temp.v_ex_m[0] = dat.le_R_ext1_m[2]
                            temp.v_in_m[0] = temp.v_ex_m[0]
                        
                        # Checar o raio do bordo de ataque (assimétricos)
                        if temp.symm_m == 0 and temp.v_ex_m[0] < dat.le_R_ext1_m[2]:
                            temp.v_ex_m[0] = dat.le_R_ext1_m[2]
                        if temp.symm_m == 0 and temp.v_in_m[0] < dat.le_R_ext2_m[2]:
                            temp.v_in_m[0] = dat.le_R_ext2_m[2]
                        
                    elif n == 2: # Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_m):
                            temp.v_ex_m[a] = dat.or_v_ex_m[a] + random.random()*dat.A_ext1_m[random.randint(0,1)]
                    
                    elif n == 3: # Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_m):
                            temp.v_in_m[a] = dat.or_v_in[a] + random.random()*dat.A_ext2_m[random.randint(0,1)]
    
                    elif n == 4: # Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                        if temp.symm_m == 1: # Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                            for a in range(1,dat.BPn_m):
                                temp.v_ex_m[a] = dat.or_v_ex_m[a] + random.random()*dat.A_ext1_m[random.randint(0,1)]
                                temp.v_in_m[a] = temp.v_ex_m[a]
                            
                        else: # Perfis assimétricos ficam com pesos intermediários distintos
                            for a in range (1,dat.BPn_m):
                                temp.v_ex_m[a] = dat.or_v_ex_m[a] + random.random()*dat.A_ext1_m[random.randint(0,1)]
                            
                            for a in range (1,dat.BPn_m):
                                temp.v_in_m[a] =dat.or_v_in[a] + random.random()*dat.A_ext2_m[random.randint(0,1)]
    
                    elif n == 5: # Alterar o ângulo do bordo de fuga
                        if temp.symm_m == 1: # Se for simétrico
                           temp.v_ex_m[dat.BPn_m] = dat.or_v_ex_m[dat.BPn_m] + random.random()*dat.B_ext1_m[random.randint(0,1)]
                           temp.v_in_m[dat.BPn_m] = temp.v_ex_m[dat.BPn_m]
                               
                           # Checar separação do bordo de fuga
                           if temp.v_ex_m[dat.BPn_m] + temp.v_in_m[dat.BPn_m] < dat.B_ext2_m[0]:
                               temp.v_ex_m[dat.BPn_m] = dat.B_ext2_m[0]/2
                               temp.v_in_m[dat.BPn_m] = dat.B_ext2_m[0]/2
                           
                        else: # Se for assimétrico
                            temp.v_ex_m[dat.BPn_m] = dat.or_v_ex_m[dat.BPn_m] + random.random()*dat.B_ext1_m[random.randint(0,1)]
                            temp.v_in_m[dat.BPn_m] = (dat.B_ext2_m[0] - temp.v_ex_m[dat.BPn_m]) + random.random()*dat.B_ext2_m[1]
                            
                            # Checar separação do bordo de fuga [beta2>=L-beta1]
                            if temp.v_in_m[dat.BPn_m] < dat.B_ext2_m[0] - temp.v_ex_m[dat.BPn_m]:
                                temp.v_in_m[dat.BPn_m] = dat.B_ext2_m[0] - temp.v_ex_m[dat.BPn_m]


                    # Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    # igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_m[1:dat.BPn_m])
                    sum2 = sum(temp.v_in_m[1:dat.BPn_m])
                    if sum2 > sum1:continue
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_m,temp.v_in_m,dat,[dat.N1_m,dat.N2_m]),dat)
            
            
                chi[s] = temp.deepcopy()
            
            # Aerofólio da ponta
            if mu_op_v[7] == 1:
                check = 0
                while check == 0:
                    
                    temp = chi[s].deepcopy()
                    if temp.symm_t == 1: # Caso o perfil seja simétrico
                        n =([1,4,5])[random.randint(0,2)]
                    else: # Caso o perfil seja assimétrico
                        n = random.randint(1,5)
                    
                    if n == 1: # Alterar o raio do bordo de ataque
                        if temp.symm_t == 1:
                            P = 1
                        else:
                            P = random.randint(1,4)
                        
                        if P == 1: # Mudar ambos para o mesmo valor
                            temp.v_ex_t[0] = dat.or_v_ex_t[0] + random.random()*dat.le_R_ext1_t[random.randint(0,1)]
                            temp.v_in_t[0] = temp.v_ex_t[0]
                        elif P == 2: # Mudar ambos independentemente 
                            temp.v_ex_t[0] = dat.or_v_ex_t[0] + random.random()*dat.le_R_ext1_t[random.randint(0,1)]
                            temp.v_in_t[0] = dat.or_v_in_t[0] + random.random()*dat.le_R_ext2_t[random.randint(0,1)]
                        elif P == 3: # Mudar do extradorso
                            temp.v_ex_t[0] = dat.or_v_ex_t[0] + random.random()*dat.le_R_ext1_t[random.randint(0,1)]
                        else: # Mudar do intradorso
                            temp.v_in_t[0] = dat.or_v_in_t[0] + random.random()*dat.le_R_ext2_t[random.randint(0,1)]
                            
                        # Checar o raio do bordo de ataque (simétricos)
                        if temp.symm_t == 1 and temp.v_ex_t[0] < dat.le_R_ext1_t[2]:
                            temp.v_ex_t[0] = dat.le_R_ext1_t[2]
                            temp.v_in_t[0] = temp.v_ex_t[0]
                        
                        # Checar o raio do bordo de ataque (assimétricos)
                        if temp.symm_t == 0 and temp.v_ex_t[0] < dat.le_R_ext1_t[2]:
                            temp.v_ex_t[0] = dat.le_R_ext1_t[2]
                        if temp.symm_t == 0 and temp.v_in_t[0] < dat.le_R_ext2_t[2]:
                            temp.v_in_t[0] = dat.le_R_ext2_t[2]
                        
                    elif n == 2: # Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_t):
                            temp.v_ex_t[a] = dat.or_v_ex_t[a] + random.random()*dat.A_ext1_t[random.randint(0,1)]
                    
                    elif n == 3: # Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                        for a in range (1,dat.BPn_t):
                            temp.v_in_t[a] = dat.or_v_in[a] + random.random()*dat.A_ext2_t[random.randint(0,1)]
    
                    elif n == 4: # Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                        if temp.symm_t == 1: # Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                            for a in range(1,dat.BPn_t):
                                temp.v_ex_t[a] = dat.or_v_ex_t[a] + random.random()*dat.A_ext1_t[random.randint(0,1)]
                                temp.v_in_t[a] = temp.v_ex_t[a]
                            
                        else: # Perfis assimétricos ficam com pesos intermediários distintos
                            for a in range (1,dat.BPn_t):
                                temp.v_ex_t[a] = dat.or_v_ex_t[a] + random.random()*dat.A_ext1_t[random.randint(0,1)]
                            
                            for a in range (1,dat.BPn_t):
                                temp.v_in_t[a] =dat.or_v_in[a] + random.random()*dat.A_ext2_t[random.randint(0,1)]
    
                    elif n == 5: # Alterar o ângulo do bordo de fuga
                        if temp.symm_t == 1: # Se for simétrico
                           temp.v_ex_t[dat.BPn_t] = dat.or_v_ex_t[dat.BPn_t] + random.random()*dat.B_ext1_t[random.randint(0,1)]
                           temp.v_in_t[dat.BPn_t] = temp.v_ex_t[dat.BPn_t]
                               
                           # Checar separação do bordo de fuga
                           if temp.v_ex_t[dat.BPn_t] + temp.v_in_t[dat.BPn_t] < dat.B_ext2_t[0]:
                               temp.v_ex_t[dat.BPn_t] = dat.B_ext2_t[0]/2
                               temp.v_in_t[dat.BPn_t] = dat.B_ext2_t[0]/2
                           
                        else: # Se for assimétrico
                            # Ignorar a mutação de beta pra perfis assimétricos caso
                            # a intenção seja não alterar o ângulo de cima. Devido à
                            # forma que o ângulo do intradorso foi implementado, ele
                            # necessariamente depende do ângulo do extradorso. O que 
                            # se deve fazer depois é permitir a manipulação do ângulo
                            # do intradorso igual ao que é feito no extradorso, mas 
                            # em seguida fazendo checagens e alterando o valor se 
                            # preciso (como os requisitos das plantas das asas).
                            # No caso do não cumprimento do requisito, deve-se alterar
                            # ou apenas um dos ângulos ou ambos igualmente.
                            if sum(dat.B_ext1_r) != 0:
                                temp.v_ex_t[dat.BPn_t] = dat.or_v_ex_t[dat.BPn_t] + random.random()*dat.B_ext1_t[random.randint(0,1)]
                                temp.v_in_t[dat.BPn_t] = (dat.B_ext2_t[0] - temp.v_ex_t[dat.BPn_t]) + random.random()*dat.B_ext2_t[1]
                                
                                # Checar separação do bordo de fuga [beta2>=L-beta1]
                                if temp.v_in_t[dat.BPn_t] < dat.B_ext2_t[0] - temp.v_ex_t[dat.BPn_t]:
                                    temp.v_in_t[dat.BPn_t] = dat.B_ext2_t[0] - temp.v_ex_t[dat.BPn_t]


                    # Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    # igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_t[1:dat.BPn_t])
                    sum2 = sum(temp.v_in_t[1:dat.BPn_t])
                    if sum2 > sum1:continue
                    
                    # Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_t,temp.v_in_t,dat,[dat.N1_t,dat.N2_t]),dat)
            
            
                chi[s] = temp.deepcopy()
            
            # Torção geométrica na ponta
            if mu_op_v[8] == 1:
                chi[s].tw_t = dat.or_tw_t + random.random()*dat.tw_t_ext(random.randint(0,1))
            
            # Enflechamento total (asas trapezoidais simples)
            if mu_op_v[9] == 1 and chi[s].type == 0 and not 'Z' in dat.sweep_ext_in:
                chi[s].sweep = dat.or_sweep + random.random()*dat.sweep_ext[random.randint(0,1)]
        
            # Enflechamento da primeira seção (asas bitrapezoidais)
            if mu_op_v[10] == 1 and chi[s].type == 1 and not 'Z' in dat.sweep1_ext_in:
                chi[s].sweep1 = dat.or_sweep1 + random.random()*dat.sweep1_ext[random.randint(0,1)]
            
            # Enflechamento da segunda seção (asas bitrapezoidais)
            if mu_op_v[11] == 1 and chi[s].type == 1 and not 'Z' in dat.sweep2_ext_in:
                chi[s].sweep2 = dat.or_sweep2 + random.random()*dat.sweep2_ext[random.randint(0,1)]
            
            # Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
            if chi[s].type == 1:
                if chi[s].b1 > chi[s].b-2*dat.b1_ext_in[2] or chi[s].c_t > chi[s].c_m and dat.or_c_m != 'L' or chi[s].c_m > chi[s].c_r and dat.or_c_m != 'L':
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
        
        
end = time()
mi = np.fix((end-start)/60)
s = (end-start)%60
print('Tempo: ' + str(mi) + 'min e ' + str(s) + 's')

# Comparar aerofólios
plt.figure()
#plot_airfoil_cst_TCC2(run_cst_TCC2(pop(pos).v_ex,pop(pos).v_in,dat),1,loop),hold on
plot_airfoil_cst_TCC2(4,run_cst_TCC2_3D(pop[pos].v_ex_r,pop[pos].v_in_r,dat,[dat.N1_r,dat.N2_r]))
plot_airfoil_cst_TCC2(3,run_cst_TCC2_3D(dat.or_v_ex_r,dat.or_v_in_r,dat,[dat.N1_r,dat.N2_r]))
plt.legend('Otimizado','Original');plt.title('Raiz');plt.axis('equal');plt.grid('True')
if dat.type == 1:
    plt.figure()
    plot_airfoil_cst_TCC2(4,run_cst_TCC2_3D(pop[pos].v_ex_m,pop[pos].v_in_m,dat,[dat.N1_m,dat.N2_m]))
    plot_airfoil_cst_TCC2(3,run_cst_TCC2_3D(dat.or_v_ex_m,dat.or_v_in_m,dat,[dat.N1_m,dat.N2_m]))
    plt.legend('Otimizado','Original');plt.title('Meio');plt.axis('equal');plt.grid('True')

plt.figure()
plot_airfoil_cst_TCC2(4,run_cst_TCC2_3D(pop[pos].v_ex_t,pop[pos].v_in_t,dat,[dat.N1_t,dat.N2_t]))
plot_airfoil_cst_TCC2(3,run_cst_TCC2_3D(dat.or_v_ex_t,dat.or_v_in_t,dat,[dat.N1_r,dat.N2_t]))
plt.legend('Otimizado','Original');plt.title('Ponta');plt.axis('equal');plt.grid('True')

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


# Pegar o struct de arquivo e imprimir todos
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
    
    print('Raiz (Polinômio de grau ' + str(dat.BPn_r) + '):')
    print('v_ex_r = [{:.4f}, '.format(archive[i].v_ex_r[0]),end='')
    for j in range(1,len(pop[0].v_ex_r)-3):
        print('{:.4f}, '.format(archive[i].v_ex_r[j]),end='')
    print('{:.4f}, '.format(archive[i].v_ex_r[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_ex_r[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_ex_r[-1]))

    print('v_in_r = [{:.4f}, '.format(archive[i].v_in_r[0]),end='')
    for j in range(1,len(pop[0].v_ex_r)-3):
        print('{:.4f}, '.format(archive[i].v_in_r[j]),end='')
    print('{:.4f}, '.format(archive[i].v_in_r[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_in_r[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_in_r[-1]))
    
    if archive[i].type == 1:
        print('Meio (Polinômio de grau ' + str(dat.BPn_m) + '):')
        print('v_ex_m = [{:.4f}, '.format(archive[i].v_ex_m[0]),end='')
        for j in range(1,len(pop[0].v_ex_m)-3):
            print('{:.4f}, '.format(archive[i].v_ex_m[j]),end='')
        print('{:.4f}, '.format(archive[i].v_ex_m[-3]),end='')
        print('{:.4f}, '.format(archive[i].v_ex_m[-2]),end='')
        print('{:.4f}];'.format(archive[i].v_ex_m[-1]))
    
        print('v_in_m = [{:.4f}, '.format(archive[i].v_in_m[0]),end='')
        for j in range(1,len(pop[0].v_ex_m)-3):
            print('{:.4f}, '.format(archive[i].v_in_m[j]),end='')
        print('{:.4f}, '.format(archive[i].v_in_m[-3]),end='')
        print('{:.4f}, '.format(archive[i].v_in_m[-2]),end='')
        print('{:.4f}];'.format(archive[i].v_in_m[-1]))
    
    print('Ponta (Polinômio de grau ' + str(dat.BPn_t) + '):')
    print('v_ex_t = [{:.4f}, '.format(archive[i].v_ex_t[0]),end='')
    for j in range(1,len(pop[0].v_ex_t)-3):
        print('{:.4f}, '.format(archive[i].v_ex_t[j]),end='')
    print('{:.4f}, '.format(archive[i].v_ex_t[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_ex_t[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_ex_t[-1]))

    print('v_in_t = [{:.4f}, '.format(archive[i].v_in_t[0]),end='')
    for j in range(1,len(pop[0].v_ex_t)-3):
        print('{:.4f}, '.format(archive[i].v_in_t[j]),end='')
    print('{:.4f}, '.format(archive[i].v_in_t[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_in_t[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_in_t[-1]))
    
    print('Torção geométrica na ponta tw_t = ' + str(archive[i].tw_t))
    
    print('- Dados aerodinâmicos -')
    for j in range(dat.cases):
#        fprintf('Condição de voo %d (CL,CD,L/D,CM): ',j),disp(archive(i).aero(j,:))
        print('Condição de voo %d: (CL,CD,L/D,CM)',j)
        print('Original: %f, %f, %f, %f'%(archive[i].aero[j,0],archive[i].aero[j,1],archive[i].aero[j,2],archive[i].aero[j,3]),end='')
        print('')
        print('Otimizado: %f, %f, %f, %f'%(original.aero[j,0],original.aero[j,1],original.aero[j,2],original.aero[j,3]))
        
        if 'q' in dat.coeff_op[:,0] or '#' in dat.coeff_op[:,0]:
            print('L = %f N (Original)'%(original.aero[j,0]*1/2*dat.rho[j]*dat.v_ref[j]**2*original.S))
            print('L = %f N (Otimizado)'%(archive[i].aero[j,0]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))
        
        if 'q' in dat.coeff_op[:,1] or '#' in dat.coeff_op[:,1]:
            print('D = %f N (Original)'%(original.aero[j,1]*1/2*dat.rho[j]*dat.v_ref[j]**2*original.S))
            print('D = %f N (Otimizado)'%(archive[i].aero[j,1]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))
        
        if 'q' in dat.coeff_op[:,3]:
            print('M = %f Nm (Original)'%(original.aero[j,3]*1/2*dat.rho[j]*dat.v_ref[j]**2*original.S))
            print('M = %f Nm (Otimizado)'%(archive[i].aero[j,3]*1/2*dat.rho[j]*dat.v_ref[j]**2*archive[i].S))
    
    print('Pontuação: %f'%(archive[i].score))
    print('')
    print('')
        
        
        