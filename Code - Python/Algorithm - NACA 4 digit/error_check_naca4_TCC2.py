from warnings import warn
from time import sleep
import numpy as np

def error_check_naca4_TCC2(dat):

    if dat.N%2 != 0:
        warn('Tamanho da população deve ser par. O valor inserido será alterado.'),sleep(5)
        dat.N += 1

    if dat.mu < 0 or dat.mu > 1:
        raise TypeError('Chance de mutação deve estar no intervalo 0 <= dat.mu <= 1')
    
    if dat.iter < 1 or dat.iter != np.floor(dat.iter):
        raise TypeError('Número de iterações do algoritmo deve ser positivo, não nulo e inteiro')

    if dat.m_ext[0] < 0 or dat.m_ext[1] > 9:
        raise TypeError('Extensão de m deve ser especificada de 0 a 9')

    if dat.p_ext[0] < 1 or dat.p_ext[1] > 9:
        raise TypeError('Extensão de p deve ser especificada de 1 a 9')

    if dat.t_ext[0] < 1 or dat.t_ext[1] > 99:
        raise TypeError('Extensão de t deve ser especificada de 1 a 99')
        
    if dat.cases < 1 or dat.cases != np.floor(dat.caes):
        raise TypeError('Número de condições de voo deve ser positivo, não nulo e inteiro')
        
    if len(dat.reynolds) < dat.cases:
        raise TypeError('Não há números de Reynolds suficientes para todas as condições de voo')
        
    if len(dat.aoa) < dat.cases:
        raise TypeError('Não há ângulos de ataque suficientes para todas as condições de voo')
        
    if len(dat.iter_sim) < dat.cases:
        raise TypeError('Não há números de iterações (XFOIL) suficientes para todas as condições de voo')
        
    # Trocar valores nulos e negativos de números de iteração pelo valor padrão do XFOIL
    for i in range(dat.cases):
        if dat.iter_sim[i] == 0:
            dat.iter_sim = 10
    
    if dat.coeff_op.shape[0] < dat.cases or dat.coeff_val.shape[0] < dat.cases or dat.coeff_F.shape[0] < dat.cases:
        raise TypeError('Configurações das funções objetivo são insuficientes para todas as condições de voo')
    
    if sum(sum(dat.coeff_op[:,0:2] == 'c')) != 0 or sum(sum(dat.coeff_op[:,0:2] == 'k')) != 0:
        raise TypeError('funções objetivo c e k valem apenas para o coeficiente de momento')

    if dat.cases == 1 and dat.coeff_op[0,3] == 'c' or dat.cases == 1 and dat.coeff_op[0,3] == 'k':
        warn('função objetivo c/k vale apenas para múltiplas condições de voo. A mesma será ignorada.')
        dat.coeff_op[0,3] = '!', sleep(5)
        
    if sum(sum(dat.coeff_op[0:dat.cases,:] == '!')) == dat.coeff_op[0:dat.cases].size:
        raise TypeError('Ao menos uma das funções objetivo deve estar ativa')
    
    for i in range(dat.cases):
        for j in range(4):
            if dat.coeff_op[i,j] != '!' and dat.coeff_op[i,j] != '^' and dat.coeff_op[i,j] != 'o' and dat.coeff_op[i,j] != 'c' and dat.coeff_op[i,j] != 'k'    :
                raise TypeError('A matriz dat.coeff_op deve ser definida com opções !, ^, o, c ou k')

    if dat.cases > 1:
        for i in range(1,dat.cases):
            if sum(dat.coeff_op[i,:] == '!') == 4 and dat.coeff_op[0,3] != 'c' and dat.coeff_op[0,3] != 'k':
                raise TypeError('Condição de voo ' + str(i+1) + ' não tem função objetivo definida')
        
    for P in range(dat.cases):
        if dat.coeff_op[P,1] == 'o' and dat.coeff_val[P,1] < 0:
            raise TypeError('Condição de voo ' + str(i) + ': CD alvo deve ser maior que zero')
        elif dat.coeff_op[P,1] == 'o' and dat.coeff_val[P,1] == 0:
            dat.coeff_op[1] = '^'
            warn('Condição de voo ' + str(i+1) + ': CD = 0 - função objetivo de arrastop trocada de o para ^'),sleep(5) 

    if dat.cases > 1 and dat.coeff_op[0,3] == 'c' or dat.cases > 1 and dat.coeff_op[0,3] == 'k':
        if sum(dat.coeff_op[1:3] != '!') != 0:
            warn('função objetivo de momento c/k: funções objetivo das condições de voo subsequentes serão ignoradas'),sleep(5)
            # temp = np.zeros((dat.cases-1,1))
            # for i in range(len(temp)):
            #     temp[i] = '!'
            dat.coeff_op[1:dat.cases,3] = np.repeat('!',dat.cases-1)
    
    
    return dat
