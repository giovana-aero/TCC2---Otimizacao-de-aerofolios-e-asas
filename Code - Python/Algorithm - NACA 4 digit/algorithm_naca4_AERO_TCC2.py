# Algoritmo genético - Perfis NACA 4 dígitos

# Lembrete: indentação automática: shift+alt+f no vscode
# Lembrete: comentários são feitos com ctrl + 1 no spyder


# a continuar depois: checar o run_xfoil_naca4_TCC2, checar o script de checagem de erros



# Pacotes necessários
import numpy as np
from ypstruct import structure
import random
# import os
import matplotlib.pyplot as plt
from time import time
import locale
locale.setlocale(locale.LC_ALL, 'pt_BR.UTF-8') # Trocar separadores decimais
plt.ticklabel_format( axis='both', style='', scilimits=None, useOffset=None, useLocale=True, useMathText=None)

# Funções utilizadas
from error_check_naca4_TCC2 import error_check_naca4_TCC2
from run_xfoil_naca4_TCC2 import run_xfoil_naca4_TCC2
from fitness_naca4 import fitness_naca4
# from make_vector import make_vector
from plot_airfoil_naca4_TCC2 import plot_airfoil_naca4_TCC2
from selection_crossover import selection_crossover

# clear = lambda: os.system('cls')
start = time() # Começar a contar o tempo de execução

# Parâmetros do algoritmo
dat = structure()                  # Inicializar o struct
dat.N = 10                             # Número de indivíduos na população
dat.mu = 0.05                           # Probabilidade de mutação
dat.iter = 3                           # Número de iterações
dat.aero = np.zeros((dat.iter,4))          # Matriz usada no gráfico final
dat.aero_m = np.zeros((dat.iter,4))
dat.elite = 1                         # Aplicar elitismo?
dat.subs = 1                          # Substituir aerofólios sem resultados?

# Parâmetros da geometria
dat.m_ext = [0,9]
dat.p_ext = [1,9]
dat.t_ext = [10,30]

# Parâmetros das simulações
dat.cases = 2
dat.reynolds = [1e6,1e6,1e6]             # Valores dos números de Reynolds para as simulações
dat.aoa = [0,2,4]                       # Ângulos de ataque
dat.iter_sim = [50,50,50]              # Números de iterações no XFOIL
dat.coeff_op = np.array([['!','o','!','k'],
                         ['!','!','!','!'],
                         ['!','!','!','!']])
dat.coeff_val = np.array([[0.6, 0.009, 90, 0],
                          [0,0,0,0],
                          [0,0,0,0]])
dat.coeff_F = np.array([[1, 1, 1, 1],
                        [1,1,1,1],
                        [1,1,1,1]])
# [CL CD L/D CM] Definição da matriz dat.coeff_op
# '!' -> não usar como função objetiva
# '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD)
# 'c'  -> buscar valor constante de coeficiente de momento (arbitrário)
# 'ce' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val[0,3])
# 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
# em sua respectiva casa na matriz dat.coeff_val
# O vetor dat.coeff_F dá os pesos de cada função objetiva

# Checagem de erros
dat = error_check_naca4_TCC2(dat)

# Gerar população inicial
empty = structure()
empty.m = None
empty.p = None
empty.t = None
empty.aero = None  # Terá o mesmo formato que o vetor coeff_op
empty.score = 0

pop = empty.repeat(dat.N)
chi = empty.repeat(dat.N)

# Gerar numeração aleatória
print('<< Geração da população inicial >>')
for i in range(dat.N):
    print('Indivíduo {}'.format(i+1))
    pop[i].m = random.randint(dat.m_ext[0], dat.m_ext[1]);  # Curvatura máxima
    if pop[i].m == 0: # Perfis simétricos
        pop[i].p = 0
    else: # Perfis assimétricos
        pop[i].p = random.randint(dat.p_ext[0], dat.p_ext[1])

    pop[i].t = random.randint(dat.t_ext[0], dat.t_ext[1])  # Espessura máxima

# Gerar struct que guarda o melhor perfil de cada geração
archive = empty.repeat(dat.iter)

# Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for loop in range(dat.iter):
    # Simular os perfis e obter dados
    select = np.ones((1,dat.N))
    print('<< Simulação dos aerofólios >>')
    for i in range(dat.N):
        print('Indivíduo {}'.format(i+1))
        naca_num = str(pop[i].m) + str(pop[i].p) + str(pop[i].t)
        # print(naca_num)
        pop[i].aero = run_xfoil_naca4_TCC2(naca_num,dat) # Simulação
        
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
    pop = fitness_naca4(pop,dat,select2)
    # # CL (aero[0])
    # if dat.coeff_op[0] != '!':
    #     F1 = dat.coeff_F[0]
    #     if dat.coeff_op[0] == '^': # Maior CL possível
    #         temp = make_vector(pop,0,select2) # Vetor com todos os CLs
    #         for i in select2:
    #             pop[i].score = pop[i].aero[0]/max(temp)*F1
    #     else: # Alcançar um CL específico
    #         CLtgt = dat.coeff_val[0]
    #         if CLtgt == 0: # Se o CL alvo for igual a zero
    #             temp = make_vector(pop,0,select2) # Vetor com todos os CLs
    #             for i in select2:
    #                 pop[i].score = (1-abs(pop[i].aero[0])/max(abs(temp)))*F1
    #         elif CLtgt > 0: # Se o CL alvo for positivo
    #             for i in select2:
    #                 if pop[i].aero[0] >= CLtgt: # Se o CL for maior/igual que o CL alvo
    #                     pop[i].score = CLtgt/pop[i].aero[0]*F1
    #                 elif pop[i].aero[0] < CLtgt: # Se o CL for menor que o CL alvo
    #                     pop[i].score = pop[i].aero[0]/CLtgt*F1
    #         else: # Se o CL alvo for negativo
    #             for i in select2:
    #                 if pop[i].aero[0] <= CLtgt: # Se o CL for menor/igual que o CL alvo
    #                     pop[i].score = CLtgt/pop[i].aero[0]*F1
    #                 else: # Se o CL for maior que o CL alvo
    #                     pop[i].score = pop[i].aero[0]/CLtgt*F1
    # # CD (aero[1])
    # if dat.coeff_op[1] != '!':
    #     F2 = dat.coeff_F[1]
    #     if dat.coeff_op[1] == '^': # Menor CD possível (tende a zero, mas nunca igual a ou menor que zero)
    #         temp = make_vector(pop,1,select2)
    #         for i in select2:
    #             if pop[i].aero[1] <= 0:
    #                 pop[i].score = 0 # Punir aerofólios com resultados irreais
    #             else:
    #                 pop[i].score = (1-pop[i].aero[1]/max(temp))*F2 + pop[i].score
    #     else: # Alcançar um CD específico
    #         CDtgt = dat.coeff_val[1]
    #         for i in select2:
    #             if pop[i].aero[1] >= CDtgt: # Se o CD for maior/igual que o CD alvo
    #                 pop[i].score = CDtgt/pop[i].aero[1]*F2 + pop[i].score
    #             else: # Se o CD for menor que o CD alvo
    #                 pop[i].score = pop[i].aero[1]/CDtgt*F2 + pop[i].score
    # # L/D (aero[2])
    # if dat.coeff_op[2] != '!':
    #     F3 = dat.coeff_F[2]
    #     if dat.coeff_op[2] == '^': # Maior L/D possível
    #         temp = make_vector(pop,2,select2) # Vetor com todos os L/Ds
    #         for i in select2:
    #             pop[i].score = pop[i].aero[2]/max(temp)*F3 + pop[i].score
    #     else: # Alcançar um L/D específico
    #         LDtgt = dat.coeff_val[2]
    #         if LDtgt == 0: # Se o L/D alvo for igual a zero
    #             temp = make_vector(pop,2,select2) # Vetor com todos os L/Ds
    #             for i in select2:
    #                 pop[i].score = (1-abs(pop[i].aero[2])/max(abs(temp)))*F3 + pop[i].score
    #         elif LDtgt > 0: # Se o L/D alvo for positivo
    #             for i in select2:
    #                 if pop[i].aero[2] >= LDtgt: # Se o L/D for maior/igual que o L/D alvo
    #                     pop[i].score = LDtgt/pop[i].aero[2]*F3 + pop[i].score
    #                 else:# Se o L/D for menor que o L/D alvo
    #                     pop[i].score = pop[i].aero[2]/LDtgt*F3 + pop[i].score
    #         else: # Se o L/D alvo for negativo
    #             if pop[i].aero[2] <= LDtgt: # Se o L/D for menor/igual que o L/D alvo
    #                 pop[i].score = LDtgt/pop[i].aero[2]*F3 + pop[i].score
    #             else: # Se o L/D for maior que o L/D alvo
    #                 pop[i].score = pop[i].aero[2]/LDtgt*F3 + pop[i].score
    # # CM (aero[3])
    # if dat.coeff_op[3] != '!':
    #     F4 = dat.coeff_F[3]
    #     if dat.coeff_op[3] == '^': # Achar o menor valor possível
    #         temp = make_vector(pop,3,select2)
    #         for i in select2:
    #             pop[i].score = (1-pop[i].aero[3]/max(temp))*F4 + pop[i].score
    #     else: # Alcançar CM específico
    #         CMtgt = dat.coeff_val[3]
    #         if CMtgt == 0: # Se o CM alvo for igual a zero
    #             temp = make_vector(pop,3,select2) # Vetor com todos os CMs
    #             for i in select2:
    #                 pop[i].score = (1-abs(pop[i].aero[3])/max(abs(temp)))*F4 + pop[i].score
    #         elif CMtgt > 0: # Se o CM alvo for positivo
    #             for i in select2:
    #                 if pop[i].aero[3] >= CMtgt: # Se o L/D for maior/igual que o L/D alvo
    #                     pop[i].score = CMtgt/pop[i].aero[3]*F4 + pop[i].score
    #                 elif pop[i].aero[3] < CMtgt: # Se o L/D for menor que o L/D alvo
    #                     pop[i].score = pop[i].aero[3]/CMtgt*F4 + pop[i].score
    #         else: # Se o CM alvo for negativo
    #             for i in select2:
    #                 if pop[i].aero[3] <= CMtgt: # Se o CM for menor/igual que o CM alvo
    #                     pop[i].score = CMtgt/pop[i].aero[3]*F4 + pop[i].score
    #                 else: # Se o CM for maior que o CM alvo
    #                     pop[i].score = pop[i].aero[3]/CMtgt*F4 + pop[i].score

    # Pôr todas as pontuações em um vetor
    weights = np.zeros(dat.N)
    for i in range(dat.N):
        weights[i] = pop[i].score
    
    # Guardar o melhor perfil de cada iteração
    pos = np.argmax(weights)
    archive[loop] = pop[pos].deepcopy()

    # Mostrar o melhor perfil
    plt.figure()
    plot_airfoil_naca4_TCC2(pop[pos],loop+1)

    # Parar o código aqui na última iteração, já que nesse cenário o resto do código é inútil
    if loop == dat.iter-1:
        break
    
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
                for i in select[0,:]:
                    pop[i] = pop[ind[k]].deepcopy()
                    k += 1
                    if k == len(ind): k = 0

    # Escolher membros da população pra reprodução
    c = 0
    print('<< Reprodução >>')
    for f in range(int(dat.N/2)):
        print('{:.2f}% completo'.format((f+1)/(dat.N/2)*100))
		#fprintf('%.2f%% completo\n',f/(dat.N/2)*100)
		
        # Isto seleciona dois pais por meio de uma seleção via roleta
        # (indivíduos com pesos maiores têm mais chance de serem selecionados)
        par = [0,0] # Vetor que indica a numeração dos pais escolhidos
        par[0] = selection_crossover(weights)
        par[1] = selection_crossover(weights)
        
        # Crossover, gerando dois filhos de cada dois pais
        if pop[par[0]].m == 0 or pop[par[1]].m == 0: # Caso um deles seja simétrico
            n = 1
        else:
            n = random.randint(1,3)
        
        if n == 1:
            chi[c].m = pop[par[0]].m
            chi[c].p = pop[par[0]].p
            chi[c].t = pop[par[1]].t
            chi[c+1].m = pop[par[1]].m
            chi[c+1].p = pop[par[1]].p
            chi[c+1].t = pop[par[0]].t

        elif n == 2:
            chi[c].m = pop[par[0]].m
            chi[c].p = pop[par[1]].p
            chi[c].t = pop[par[0]].t
            chi[c+1].m = pop[par[1]].m
            chi[c+1].p = pop[par[0]].p
            chi[c+1].t = pop[par[1]].t
                
        else:
            chi[c].m = pop[par[1]].m
            chi[c].p = pop[par[0]].p
            chi[c].t = pop[par[0]].t
            chi[c+1].m = pop[par[0]].m
            chi[c+1].p = pop[par[1]].p
            chi[c+1].t = pop[par[1]].t           

        c += 2
    
    # Mutação
    select1 = int(np.random.uniform(0,1,len(pop)) <= dat.mu)
    #select1 = (rand(size(pop)) <= dat.mu)
    if sum(select1) != 0:
        print('<< Mutação >>')
        select2 = np.argwhere(select1 == True)
        k = 0
        for i in select2:
            s = select2[k,0]
            
            
            # Código da mutação referente apenas aos perfis NACA 4 dígitos
            # começa aqui
            if chi[int(s)].m == 0: # Perfis simétricos
                # random.randint(dat.p_ext[0], dat.p_ext[1])
                n = random.randint(1,2)
                if n == 1: # Transformar em um perfil assimétrico
                    chi[s].m = random.randint(1,dat.m_ext[1])
                    chi[s].p = random.randint(dat.p_ext[0],dat.p_ext[1])
                else: # Mudar apenas a espessura
                    chi[s].t = random.randint(dat.t_ext[0],dat.t_ext[1])

            else: # Perfis assimétricos
                n = random.randint(1,4)
                if n == 1: # Mudar apenas a curvatura máxima
                    chi[s].m = random.randint(1,dat.m_ext[1])
                elif n == 2: # Mudar apenas o local da curvatura máxima
                    chi[s].p = random.randint(dat.p_ext[0],dat.p_ext[1])
                elif n == 3: # Mudar apenas a espessura
                    chi[s].t = random.randint(dat.t_ext[0],dat.t_ext[1])
                else: # Transformar em um perfil simétrico
                    chi[s].m = 0
                    chi[s].p = 0
                
            # Código da mutação referente apenas aos perfis NACA 4 dígitos
            # termina aqui
            
            k += 1

    # Substituir a população inicial pelos filhos
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
    
# Pegar o struct de arquivo de imprimir todos
for i in range(len(archive)):
    print('<< Iteração {} >>'.format(i+1))
    print('NACA %d%d%d'%(archive[i].m,archive[i].p,archive[i].t))
    print('- Dados aerodinâmicos -')
    for j in range(dat.cases):
        print('Condição de voo %d: %f, %f, %f, %f'%(j,archive[i].aero[j,0],archive[i].aero[j,1],archive[i].aero[j,2],archive[i].aero[j,3]))
    print('Pontuação: {}'.format(archive[i].score))
    print('')
