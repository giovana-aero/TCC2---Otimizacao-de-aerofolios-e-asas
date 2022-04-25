## Algoritmo genético CST - Otimizador de perfis específicos

# Pacotes necessários
import numpy as np
from ypstruct import structure
import random
# import os
import matplotlib.pyplot as plt
from time import time
# Funções utilizadas
from error_check_cst_TCC2 import error_check_cst_TCC2
from quality import quality
from run_cst_TCC2 import run_cst_TCC2
from run_xfoil_cst_TCC2 import run_xfoil_cst_TCC2
from fitness_cst import fitness_cst
# from make_vector import make_vector
from selection_crossover import selection_crossover
from plot_airfoil_cst_TCC2 import plot_airfoil_cst_TCC2

# Começar a contar o tempo
start = time()

# Parâmetros do algoritmo
dat = structure()
dat.N = 300                              # Número de indivíduos na população
dat.mu = 0.05                           # Probabilidade de mutação
dat.iter = 3                           # Número de iterações
dat.elite = 1                         # Aplicar elitismo?
dat.subs = 1                          # Substituir aerofólios sem resultados?

# Dados do aerofólio a ser otimizado 
dat.or_v_ex = [0.020000, 0.246667, 0.253333, 0.212000, 0.158924, 0.174621, 21.000000, 0.000000]
dat.or_v_in = [0.020000, 0.066667, 0.021333, 0.022000, 0.034667, 0.033333, -0.000000, -0.000000]
dat.symm_override = 0 # Forçar a identificação de perfil simétrico caso o algoritmo não o reconheça como tal (devido a erros de precisão)
                      # Isto efetivamente forma um perfil simétrico fazendo o espelhamento do extradorso

# Parâmetros da geometria CST
dat.BPn = len(dat.or_v_ex)-2                  # Grau do polinômio de Bernstein (número de variáveis de design = BPn+1, desconsiderando o delta_z)
dat.np = 80                          # Número de pontos a serem usados na geração de ordenadas
dat.p_op = 0                # Opção da geração de pontos (1 pra cosspace, outro valor pra cosspace_half)
dat.N1 = 0.5
dat.N2 = 1
dat.le_R_ext1 = [-0.02,0.02,0.005] # Limite inferior (<=0), limite superior(>=0), valor mínimo possível (pois Rle > 0)
dat.le_R_ext2 = [-0.02,0.02,0.005] # Limite inferior (<=0), limite superior(>=0), valor mínimo possível (pois Rle > 0)
dat.A_ext1 = [-0.05,0.05] # Limite inferior (<=0) e superior (>=0)
dat.A_ext2 = [-0.05,0.05] # Limite inferior (<=0) e superior (>=0)
dat.B_ext1 = [-5,5] # Limite inferior (<=0) e superior (>=0)
dat.B_ext2 = [dat.or_v_ex[-2] + dat.or_v_in[-2],10] # O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
dat.chord = 1

# Parâmetros das simulações
dat.cases = 3                           # Número de condições de voo a serem analisadas
dat.reynolds = [1e5,1e5,1e5]            # Valores do números de Reynolds para as simulações
dat.aoa = [0,2,4]                       # Ângulos de ataque
dat.iter_sim = [20,20,10]                       # Número de iterações no XFOIL
dat.numNodes = 0                         # Número de painéis
dat.coeff_op = np.array([['!','!','!','c'],
                         ['!','!','!','o'],
                         ['!','!','!','o']])
dat.coeff_val = np.array([[0.35,0.007,90,-0.1],
                          [0.4,0,0,0],
                          [0,0,0,0]])
dat.coeff_F = np.array([[1,1,1,1],
                        [1,1,1,1],
                        [1,1,1,1]])
# [CL CD L/D CM] Definição do vetor dat.coeff_op
# '!' -> não usar como função objetiva
# '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD)
# 'c' -> buscar valor constante de coeficiente de momento (arbitrário)
# 'k' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val(1,4))
# 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
# em sua respectiva casa no vetor dat.coeff_val
# A matriz dat.coeff_F dá os pesos de cada função objetiva 

# Checagem de erros
dat = error_check_cst_TCC2(dat,2)
                                
# Template dos structs
empty = structure()
empty.v_ex = np.zeros(dat.BPn+2)
empty.v_in = np.zeros(dat.BPn+2)
empty.symm = None
empty.aero = np.zeros([dat.cases,4])
empty.score = 0

pop = empty.repeat(dat.N)
chi = empty.repeat(dat.N)

# Vetor que define os perfis:
# v = [ RLe A1 A2 A3 ... A(N) beta Dz]

# Gerar população inicial
print('<< Geração da população inicial >>')
for i in range(dat.N):
    print('Indivíduo {}'.format(i+1))
    
    check = 0
    while check == 0:
    
        if dat.or_v_ex == dat.or_v_in or dat.symm_override == 1: # Perfil simétrico
            # Vetor com informações do extradorso
            pop[i].v_ex[0] = dat.or_v_ex[0] + random.random()*dat.le_R_ext1[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn):
                pop[i].v_ex[a] = dat.or_v_ex[a] + random.random()*dat.A_ext1[random.randint(0,1)] # Pesos intermediários
            pop[i].v_ex[dat.BPn] = dat.or_v_ex[dat.BPn] + random.random()*dat.B_ext1[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex[dat.BPn+1] = dat.or_v_ex[dat.BPn+1] #randi[dat.delta_range]*0.1*rand;   # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in = pop[i].v_ex
        
            # Checar o raio do bordo de ataque
            if pop[i].v_ex[0] < dat.le_R_ext1[2]:
                pop[i].v_ex[0] = dat.le_R_ext1[2]
                pop[i].v_in[0] = pop[i].v_ex[0]
            
            # Checar separação do bordo de fuga
            if pop[i].v_ex[dat.BPn] + pop[i].v_in[dat.BPn] < dat.B_ext2[0]:
                pop[i].v_ex[dat.BPn] = dat.B_ext2[0]/2
                pop[i].v_in[dat.BPn] = dat.B_ext2[0]/2
            
            pop[i].symm = 1
            
        else: # perfil assimétrico
            # Vetor com informações do extradorso
            pop[i].v_ex[0] = dat.or_v_ex[0] + random.random()*dat.le_R_ext1[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn):
                pop[i].v_ex[a] = dat.or_v_ex[a] + random.random()*dat.A_ext1[random.randint(0,1)] # Pesos intermediários
            pop[i].v_ex[dat.BPn] = dat.or_v_ex[dat.BPn] + random.random()*dat.B_ext1[random.randint(0,1)]     # Ângulo do bordo de fuga
            pop[i].v_ex[dat.BPn+1] = dat.or_v_ex[dat.BPn+1] #randi[dat.delta_range]*0.1*rand;   # delta_z
            
            # Vetor com informações do intradorso
            pop[i].v_in[0] = dat.or_v_in[0] + random.random()*dat.le_R_ext2[random.randint(0,1)]    # Raio do bordo de ataque
            for a in range(1,dat.BPn):
                pop[i].v_in[a] = dat.or_v_in[a] + random.random()*dat.A_ext2[random.randint(0,1)] # Pesos intermediários
            pop[i].v_in[dat.BPn] = (dat.B_ext2[0] - pop[i].v_ex[dat.BPn]) + random.random()*dat.B_ext2[1]     # Ângulo do bordo de fuga
            pop[i].v_in[dat.BPn+1] = dat.or_v_in[dat.BPn+1] #randi[dat.delta_range]*0.1*rand;   # delta_z
            
            # Checar o raio do bordo de ataque
            if pop[i].v_ex[0] < dat.le_R_ext1[2]:
                pop[i].v_ex[0] = dat.le_R_ext1[2]
            
            if pop[i].v_in[0] < dat.le_R_ext2[2]:
                pop[i].v_in[0] = dat.le_R_ext2[2]
            
            # Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop[i].v_in[dat.BPn] < dat.B_ext2[0] - pop[i].v_ex[dat.BPn]:
                pop[i].v_in[dat.BPn] = dat.B_ext2[0] - pop[i].v_ex[dat.BPn]
            
            # Checar os pesos (soma de pesos do intradorso deve ser menor ou
            # igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop[i].v_ex[1:dat.BPn])
            sum2 = sum(pop[i].v_in[1:dat.BPn])
            if sum2 > sum1:continue
            
            pop[i].symm = 0
        
        
        # Checagem de qualidade
        check = quality(run_cst_TCC2(pop[i].v_ex,pop[i].v_in,dat),dat)
        
# Gerar struct que guarda o melhor perfil de cada geração
archive = empty.repeat(dat.iter)

# Montar arquivo de input pro XFOIL
run_xfoil_cst_TCC2(dat,1)

## Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for loop in range(dat.iter):
    
    # Simular os perfis e obter dados
    select = np.ones((1,dat.N))
    print('<< Simulação dos aerofólios >>')
    for i in range(dat.N):
        print('Indivíduo {}'.format(i+1))
        run_cst_TCC2(pop[i].v_ex,pop[i].v_in,dat,1)
        pop[i].aero = run_xfoil_cst_TCC2(dat,2) # Simulação
        
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
        print('Nenhum aerofólio convergiu nas simulações')
    
    # Atribuir pontuações (fitnesses)
    pop = fitness_cst(pop,dat,select2)

    # Pôr todas as pontuações em um vetor
    weights = np.zeros(dat.N)
    for i in range(dat.N):
        weights[i] = pop[i].score
    
    # Guardar o melhor perfil de cada iteração
    pos = np.argmax(weights)
    archive[loop] = pop[pos].deepcopy()

    # Mostrar o melhor perfil
    plt.figure()
    if dat.cases == 1: # Mostrar no título o número da iteração e os dados aerodinâmicos
        plot_airfoil_cst_TCC2(2,run_cst_TCC2(pop[pos].v_ex,pop[pos].v_in,dat),loop+1,pop[pos])
    else: # Mostrar no título apenas o número da iteração
        plot_airfoil_cst_TCC2(1,run_cst_TCC2(pop[pos].v_ex,pop[pos].v_in,dat),loop+1)
    plot_airfoil_cst_TCC2(3,run_cst_TCC2(dat.or_v_ex,dat.or_v_in,dat))
    plt.axis('equal');plt.grid('True')

    # Parar o código aqui na última iteração, já que nesse cenário o resto do código é inútil
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
    for f in range(int(dat.N/2)):
        print('{:.2f}% completo'.format((f+1)/(dat.N/2)*100))
        
        # Crossover, gerando dois filhos de cada dois pais
        # Aqui ignora-se os últimos dois parâetros (f e delta_z) porque eles
        # sempre serão iguais
        check = 0
        while check == 0:
            
            # Isto seleciona dois pais por meio de uma seleção via roleta
            # (indivíduos com pesos maiores têm  mais chance de serem selecionados)
            par = [0,0] # Vetor que indica a numeração dos pais escolhidos
            par[0] = selection_crossover(weights)
            par[1] = selection_crossover(weights)
            
            # v = [ RLe A1 A2 A3 ... A(N) beta f Dz]
            if pop[par[0]].symm == 0 and pop[par[1]].symm == 0: # Se ambos forem assimétricos
            
            
                n = random.randint(1,4)
                if n == 1: # Trocar os extradorsos e intradorsos inteiros
                    chi[c].v_ex = pop[par[0]].v_ex
                    chi[c].v_in = pop[par[1]].v_in
                    chi[c+1].v_ex = pop[par[1]].v_ex
                    chi[c+1].v_in = pop[par[0]].v_in
                    
                elif n == 2: # Trocar o raio do bordo de ataque
                    chi[c].v_ex = np.hstack((pop[par[0]].v_ex[0],pop[par[1]].v_ex[1:]))
                    chi[c].v_in = np.hstack((pop[par[0]].v_in[0],pop[par[1]].v_in[1:]))
                    chi[c+1].v_ex = np.hstack((pop[par[1]].v_ex[0],pop[par[0]].v_ex[1:]))
                    chi[c+1].v_in = np.hstack((pop[par[1]].v_in[0],pop[par[0]].v_in[1:]))
    
                elif n == 3: # Trocar os pesos intermediários
                    if dat.BPn == 2:
                            op = 1
                    else:
                        op = random.randint(1,2)
                    
                    if op == 1: # Trocar tudo    
                        chi[c].v_ex = np.hstack(([pop[par[1]].v_ex[0],pop[par[0]].v_ex[1:dat.BPn],pop[par[1]].v_ex[dat.BPn:]]))
                        chi[c].v_in = np.hstack((pop[par[1]].v_in[0],pop[par[0]].v_in[1:dat.BPn],pop[par[1]].v_in[dat.BPn:]))
                        chi[c+1].v_ex = np.hstack((pop[par[0]].v_ex[0],pop[par[1]].v_ex[1:dat.BPn],pop[par[0]].v_ex[dat.BPn:]))
                        chi[c+1].v_in = np.hstack((pop[par[0]].v_in[0],pop[par[1]].v_in[1:dat.BPn],pop[par[0]].v_in[dat.BPn:]))
    
                    else: # Trocar cortes
                        num1 = random.randint(1,dat.BPn-1)
                        num2 = random.randint(1,dat.BPn-1)
                        temp1_1 = np.hstack((pop[par[0]].v_ex[1:num1+1],pop[par[1]].v_ex[num1+1:dat.BPn]))
                        temp1_2 = np.hstack((pop[par[0]].v_in[1:num2+1],pop[par[1]].v_in[num2+1:dat.BPn]))
                        temp2_1 = np.hstack((pop[par[1]].v_ex[1:num1+1],pop[par[0]].v_ex[num1+1:dat.BPn]))
                        temp2_2 = np.hstack((pop[par[1]].v_in[1:num2+1],pop[par[0]].v_in[num2+1:dat.BPn]))
                        chi[c].v_ex = np.hstack(([pop[par[0]].v_ex[0],temp2_1,pop[par[0]].v_ex[dat.BPn:]]))
                        chi[c].v_in = np.hstack(([pop[par[0]].v_in[0],temp2_2,pop[par[0]].v_in[dat.BPn:]]))
                        chi[c+1].v_ex = np.hstack(([pop[par[1]].v_ex[0],temp1_1,pop[par[1]].v_ex[dat.BPn:]]))
                        chi[c+1].v_in = np.hstack(([pop[par[1]].v_in[0],temp1_2,pop[par[1]].v_in[dat.BPn:]]))
                    
                else: # Trocar os ângulos do bordo de fuga
                    chi[c].v_ex = np.hstack(([pop[par[1]].v_ex[0:dat.BPn],pop[par[0]].v_ex[dat.BPn],pop[par[1]].v_ex[dat.BPn+1:]]))
                    chi[c].v_in = np.hstack(([pop[par[1]].v_in[0:dat.BPn],pop[par[0]].v_in[dat.BPn],pop[par[1]].v_in[dat.BPn+1:]]))
                    chi[c+1].v_ex = np.hstack(([pop[par[0]].v_ex[0:dat.BPn],pop[par[1]].v_ex[dat.BPn],pop[par[0]].v_ex[dat.BPn+1:]]))
                    chi[c+1].v_in = np.hstack(([pop[par[0]].v_in[0:dat.BPn],pop[par[1]].v_in[dat.BPn],pop[par[0]].v_in[dat.BPn+1:]]))
                
                if n == 1:
                    # Checar a separação dos bordos de fuga. Se não cumprirem o 
                    # requisito de separação, alterar o ângulo do bordo de fuga
                    # do intradorso
                    if chi[c].v_in[dat.BPn] < [dat.B_ext2[1]-chi[c].v_ex[dat.BPn]]:
                        chi[c].v_in[dat.BPn] = dat.B_ext2[1]-chi[c].v_ex[dat.BPn]
                    
                    if chi[c+1].v_in[dat.BPn] < [dat.B_ext2[1]-chi[c+1].v_ex[dat.BPn]]:
                        chi[c+1].v_in[dat.BPn] = dat.B_ext2[1]-chi[c+1].v_ex[dat.BPn]
                    
                    # Decisão de alterar o intradorso em base de uma nota na página
                    # 57(87) do Raymer (2018)
                
                # Checar os pesos
                sum1 = sum(chi[c].v_ex[1:dat.BPn])
                sum2 = sum(chi[c].v_in[1:dat.BPn])
                if sum2 > sum1:continue
                sum1 = sum(chi[c+1].v_ex[1:dat.BPn])
                sum2 = sum(chi[c+1].v_in[1:dat.BPn])
                if sum2 > sum1:continue
                
                # Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                # for usado como variável)
                #chi(c).v_in(dat.BPn) = -chi(c).v_ex(dat.BPn)*chi(c).v_ex(dat.BPn-1)/chi(c).v_in(dat.BPn-1);
                #chi(c+1).v_in(dat.BPn) = -chi(c+1).v_ex(dat.BPn)*chi(c+1).v_ex(dat.BPn-1)/chi(c+1).v_in(dat.BPn-1);
                
                chi[c].symm = 0
                chi[c+1].symm = 0
                
            elif pop[par[0]].symm == 1 and pop[par[1]].symm == 1: # Se ambos forem simétricos
                n = random.randint(1,3)
                if n == 1: # Trocar o raio do bordo de ataque
                    chi[c].v_ex = np.hstack((pop[par[0]].v_ex[0],pop[par[1]].v_ex[1:]))
                    chi[c].v_in = chi[c].v_ex
                    chi[c+1].v_ex = np.hstack((pop[par[1]].v_ex[0],pop[par[0]].v_ex[1:]))
                    chi[c+1].v_in = chi[c+1].v_ex
                    
                if n == 2: # Trocar os pesos intermediários
                    if dat.BPn == 2:
                        op = 1
                    else:
                        op = random.randint(1,2)
                    
                    if op == 1: # Trocar tudo    
                        chi[c].v_ex = np.hstack((pop[par[1]].v_ex[0],pop[par[0]].v_ex[1:dat.BPn],pop[par[1]].v_ex[dat.BPn:]))
                        chi[c].v_in = chi[c].v_ex
                        chi[c+1].v_ex = np.hstack((pop[par[0]].v_ex[0],pop[par[1]].v_ex[1:dat.BPn],pop[par[0]].v_ex[dat.BPn:]))
                        chi[c+1].v_in = chi[c+1].v_ex
                        
                    else: # Trocar cortes
                        num1 = random.randint(2,dat.BPn) 
                        temp1_1 = np.hstack((pop[par[0]].v_ex[1:num1],pop[par[1]].v_ex[num1:dat.BPn]))
                        temp2_1 = np.hstack((pop[par[1]].v_ex[1:num1],pop[par[0]].v_ex[num1:dat.BPn]))
                        chi[c].v_ex = np.hstack((pop[par[0]].v_ex[0],temp2_1,pop[par[0]].v_ex[dat.BPn:]))
                        chi[c].v_in = chi[c].v_ex
                        chi[c+1].v_ex = np.hstack((pop[par[1]].v_ex[0],temp1_1,pop[par[1]].v_ex[dat.BPn:]))
                        chi[c+1].v_in = chi[c+1].v_ex
                    
                if n == 3: # Trocar os ângulos do bordo de fuga
                    chi[c].v_ex = np.hstack((pop[par[1]].v_ex[0:dat.BPn],pop[par[0]].v_ex[dat.BPn],pop[par[1]].v_ex[dat.BPn+1]))
                    chi[c].v_in = np.hstack((pop[par[1]].v_in[0:dat.BPn],pop[par[0]].v_in[dat.BPn],pop[par[1]].v_in[dat.BPn+1]))
                    chi[c+1].v_ex = np.hstack((pop[par[0]].v_ex[0:dat.BPn],pop[par[1]].v_ex[dat.BPn],pop[par[0]].v_ex[dat.BPn+1]))
                    chi[c+1].v_in = np.hstack((pop[par[0]].v_in[0:dat.BPn],pop[par[1]].v_in[dat.BPn],pop[par[0]].v_in[dat.BPn+1]))
                
                chi[c].symm = 1
                chi[c+1].symm = 1
            
            else: # Se um for simétrico e o outro for assimétrico
                # Transformar o simétrico em um assimétrico e vice-versa
                # chi(c) é simétrico e chi(c+1) é assimétrico
                if pop[par[0]].symm == 0:
                    chi[c].v_ex = pop[par[0]].v_ex
                    chi[c].v_in = pop[par[0]].v_ex
                    chi[c+1].v_ex = pop[par[1]].v_ex
                    chi[c+1].v_in = pop[par[0]].v_in
                else:
                    chi[c].v_ex = pop[par[1]].v_ex
                    chi[c].v_in = pop[par[1]].v_ex
                    chi[c+1].v_ex = pop[par[0]].v_ex
                    chi[c+1].v_in = pop[par[1]].v_in
                         
                # Checar a separação dos bordos de fuga. Se não cumprirem o 
                # requisito de separação, alterar o ângulo do bordo de fuga
                # do intradorso
                if chi[c].v_in[dat.BPn] < [dat.B_ext2[0]-chi[c].v_ex[dat.BPn]]:
                    chi[c].v_ex[dat.BPn] = dat.B_ext2[0]/2
                    chi[c].v_in[dat.BPn] = dat.B_ext2[0]/2
                
                if chi[c+1].v_in[dat.BPn] < [dat.B_ext2[0]-chi[c+1].v_ex[dat.BPn]]:
                    chi[c+1].v_in[dat.BPn] = dat.B_ext2[0]-chi[c+1].v_ex[dat.BPn]
                
                # Decisão de alterar o intradorso em base de uma nota na página
                # 57(87) do Raymer (2018)
                
                # Checar os pesos
                sum1 = sum(chi[c+1].v_ex[1:dat.BPn])
                sum2 = sum(chi[c+1].v_in[1:dat.BPn])
                if sum2 > sum1:continue
                
                chi[c].symm = 1
                chi[c+1].symm = 0
            
            # Checagem de qualidade
            check = quality(run_cst_TCC2(chi[c].v_ex,chi[c].v_in,dat),dat)
            if check == 0:continue
            check = quality(run_cst_TCC2(chi[c+1].v_ex,chi[c+1].v_in,dat),dat)

        c += 2
    
    # Mutação
    select1 = (np.random.uniform(0,1,len(pop)) <= dat.mu)
    if sum(select1) != 0:
        print('<< Mutação >>')
        select2 = np.argwhere(select1 == True)
        k = 0
        for i in select2:
            s = select2[k,0]
            
            
            
            check = 0
            while check == 0:
                
                temp = chi[s].deepcopy()
                if temp.symm == 1: # Caso o perfil seja simétrico
                    n = ([1,4,5])[random.randint(0,2)]
                else: # Caso o perfil seja assimétrico
                    n = random.randint(1,5)
                
                if n == 1: # Alterar o raio do bordo de ataque
                    if temp.symm == 1:
                        p = 1
                    else:
                        p = random.randint(1,4)
                        
                    if p == 1: # Mudar ambos para o mesmo valor
                        temp.v_ex[0] = dat.or_v_ex[0] + random.random()*dat.le_R_ext1[random.randint(0,1)]
                        temp.v_in[0] = temp.v_ex[0]
                    elif p == 2: # Mudar ambos independentemente
                        temp.v_ex[0] = dat.or_v_ex[0] + random.random()*dat.le_R_ext1[random.randint(0,1)]
                        temp.v_in[0] = dat.or_v_in[0] + random.random()*dat.le_R_ext2[random.randint(0,1)]
                    elif p == 3: # Mudar do extradorso
                        temp.v_ex[0] = dat.or_v_ex[0] + random.random()*dat.le_R_ext1[random.randint(0,1)]
                    else: # Mudar do intradorso
                        temp.v_in[0] = dat.or_v_in[0] + random.random()*dat.le_R_ext2[random.randint(0,1)]
                    
                    # Checar o raio do bordo de ataque (simétricos)
                    if temp.symm == 1 and temp.v_ex[0] < dat.le_R_ext1[2]:
                        temp.v_ex[0] = dat.le_R_ext1[2]
                        temp.v_in[0] = temp.v_ex[0]
                    
                    # Checar o raio do bordo de ataque (assimétricos)
                    if temp.symm == 0 and temp.v_ex[0] < dat.le_R_ext1[2]:
                        temp.v_ex[0] = dat.le_R_ext1[2]
                    if temp.symm == 0 and temp.v_in[0] < dat.le_R_ext2[2]:
                        temp.v_in[0] = dat.le_R_ext2[2]

                elif n == 2: # Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                    for a in range (1,dat.BPn):
                        temp.v_ex[a] = dat.or_v_ex[a] + random.random()*dat.A_ext1[random.randint(0,1)]
                    
                elif n == 3: # Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                    for a in range (1,dat.BPn):
                        temp.v_in[a] = dat.or_v_in[a] + random.random()*dat.A_ext2[random.randint(0,1)]

                elif n == 4: # Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                    if temp.symm == 1: # Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                        for a in range(1,dat.BPn):
                            temp.v_ex[a] = dat.or_v_ex[a] + random.random()*dat.A_ext1[random.randint(0,1)]
                            temp.v_in[a] = temp.v_ex[a]
                        
                    else: # Perfis assimétricos ficam com pesos intermediários distintos
                        for a in range (1,dat.BPn):
                            temp.v_ex[a] = dat.or_v_ex[a] + random.random()*dat.A_ext1[random.randint(0,1)]
                        
                        for a in range (1,dat.BPn):
                            temp.v_in[a] =dat.or_v_in[a] + random.random()*dat.A_ext2[random.randint(0,1)]

                elif n == 5: # Alterar o ângulo do bordo de fuga
                    if temp.symm == 1: # Se for simétrico
                       temp.v_ex[dat.BPn] = dat.or_v_ex[dat.BPn] + random.random()*dat.B_ext1[random.randint(0,1)]
                       temp.v_in[dat.BPn] = temp.v_ex[dat.BPn]
                           
                       # Checar separação do bordo de fuga
                       if temp.v_ex[dat.BPn] + temp.v_in[dat.BPn] < dat.B_ext2[0]:
                           temp.v_ex[dat.BPn] = dat.B_ext2[0]/2
                           temp.v_in[dat.BPn] = dat.B_ext2[0]/2
                       
                    else: # Se for assimétrico
                        temp.v_ex[dat.BPn] = dat.or_v_ex[dat.BPn] + random.random()*dat.B_ext1[random.randint(0,1)]
                        temp.v_in[dat.BPn] = (dat.B_ext2[0] - temp.v_ex[dat.BPn]) + random.random()*dat.B_ext2[1]
                        
                        # Checar separação do bordo de fuga [beta2>=L-beta1]
                        if temp.v_in[dat.BPn] < dat.B_ext2[0] - temp.v_ex[dat.BPn]:
                            temp.v_in[dat.BPn] = dat.B_ext2[0] - temp.v_ex[dat.BPn]

                # Checar os pesos (soma de pesos do intradorso deve ser menor ou
                # igual à soma de pesos do extradorso)
                sum1 = sum(temp.v_ex[1:dat.BPn])
                sum2 = sum(temp.v_in[1:dat.BPn])
                if sum2 > sum1:continue
                
                # Checagem de qualidade
                check = quality(run_cst_TCC2(temp.v_ex,temp.v_in,dat),dat)
            
            
            chi[s] = temp.deepcopy()
            k += 1
        

    # Substituir a população inicial pelos filhos
    for i in range(dat.N):
        pop[i] = chi[i].deepcopy()

    # Aplicar elitismo
    if dat.elite == 1:
        # Passar o melhor indivíduo pra nova população
        pop[0] = archive[loop].deepcopy()
        pop[0].score = 0


## Final ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fazer simulação do aerofólio original para fins de comparação
original = empty.deepcopy()
run_cst_TCC2(dat.or_v_ex,dat.or_v_in,dat,1)
original.aero = run_xfoil_cst_TCC2(dat,2)

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

# Pegar o struct de arquivo e imprimir todos
print('Grau do polinômio: {}'.format(dat.BPn))
for i in range(len(archive)):

    print('<< Iteração {} >>'.format(i+1))
    
    print('v_ex = [{:.4f}, '.format(archive[i].v_ex[0]),end='')
    for j in range(1,len(pop[0].v_ex)-3):
        print('{:.4f}, '.format(archive[i].v_ex[j]),end='')
    print('{:.4f}, '.format(archive[i].v_ex[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_ex[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_ex[-1]))

    print('v_in = [{:.4f}, '.format(archive[i].v_in[0]),end='')
    for j in range(1,len(pop[0].v_ex)-3):
        print('{:.4f}, '.format(archive[i].v_in[j]),end='')
    print('{:.4f}, '.format(archive[i].v_in[-3]),end='')
    print('{:.4f}, '.format(archive[i].v_in[-2]),end='')
    print('{:.4f}];'.format(archive[i].v_in[-1]))
        
    
    print('- Dados aerodinâmicos ')
    for j in range(dat.cases):
       print('Condição de voo %d: (CL,CD,L/D,CM)',j)
       print('Original: %f, %f, %f, %f'%(archive[i].aero[j,0],archive[i].aero[j,1],archive[i].aero[j,2],archive[i].aero[j,3]),end='')
       print('')
       print('Otimizado: %f, %f, %f, %f'%(original.aero[j,0],original.aero[j,1],original.aero[j,2],original.aero[j,3]))
    print('Pontuação: {}'.format(archive[i].score))
    print('')


# Template de comando pra pegar coordenadas e um gráfico de determinado perfil
#plot_airfoil_cst_TCC2(1,run_cst_TCC2(v_ex,v_in,dat,0))

#figure(2),saveas(gcf,'fig1.png')
#figure(3),saveas(gcf,'fig2.png')


