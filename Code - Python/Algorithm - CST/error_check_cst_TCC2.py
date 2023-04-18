from warnings import warn
from time import sleep
from numpy import floor,repeat

# op = 1 pro algoritmo normal
# op = 2 pro algoritmo específico

def error_check_cst_TCC2(dat,op=1):

    if dat.N%2 != 0:
        warn('Tamanho da população deve ser par. O valor inserido será alterado');sleep(5)
        dat.N += 1

    if dat.N1 != 0.5 or dat.N2 != 1:
        warn('Coeficientes N1 e N2 da função class estão especificando uma geometria que não é um aerofólio')
        sleep(5)

    if dat.mu < 0 or dat.mu > 1:
        raise TypeError('Chance de mutação deve estar no intervalo 0 <= dat.mu <= 1');sleep(5)
    
    if dat.iter < 1 or dat.iter != floor(dat.iter):
        raise TypeError('Número de iterações do algoritmo deve ser positivo, não nulo e inteiro')

    if op == 1: # Válido apenas pro algoritmo normal
        if dat.symm_op < 0 or dat.symm_op > 1:
            raise TypeError('Opção de perfis simétricos deve estar definida no intervalo 0 <= dat.symm_op <= 1')

    if dat.cases < 1 or dat.cases != floor(dat.cases):
        raise TypeError('Número de condições de voo deve ser positivo, não nulo e inteiro')

    if len(dat.reynolds) < dat.cases:
        raise TypeError('Não há números de Reynolds suficientes para todas as condições de voo')

    if len(dat.aoa) < dat.cases:
        raise TypeError('Não há ângulos de ataque suficientes para todas as condições de voo')
        
    if len(dat.iter_sim) < dat.cases:
        raise TypeError('Não há números de iterações (XFOIL) suficientes para todas as condições de voo')
        
    # Trocar valores nulos e negativos de números de iterações pelo valor padrão do XFOIL
    for i in range(dat.cases):
        if dat.iter_sim[i] <= 0:
            dat.iter_sim[i] = 10
    
    if dat.coeff_op.shape[0] < dat.cases or dat.coeff_val.shape[0] < dat.cases or dat.coeff_F.shape[0] < dat.cases:
        raise TypeError('Configurações das funções objetivo são insuficientes para todas as condições de voo')
        
    if 'c' in dat.coeff_op[:,0:3] or 'k' in dat.coeff_op[:,0:3]:
        raise TypeError('funções objetivo c e k valem apenas para o coeficiente de momento')
    
    if dat.cases == 1 and dat.coeff_op[0,3] == 'c' or dat.cases == 1 and dat.coeff_op[0,3] == 'k':
        warn('função objetivo c/k vale apenas para múltiplas condições de voo. A mesma será ignorada')
        dat.coeff_op[0,3] = '!'; sleep(5)
    
    if sum(sum(dat.coeff_op=='!')) == dat.coeff_op.size:
        raise TypeError('Ao menos uma das funções objetivo deve estar ativa')

    if op == 2: # Válida apenas pro algoritmo específico
        if dat.or_v_ex[dat.BPn] + dat.B_ext2[1] < dat.B_ext2[0]:
            raise TypeError('Extensão de Beta 2 insificiente para cumprir o requisito de separação')
    
    if op == 2: # Válida apenas pro algoritmo específico
        if dat.le_R_ext1[2] <= 0 or dat.le_R_ext2[2] <= 0:
            raise TypeError('Valor mínimo do raio do bordo de ataque deve ser maior do que zero')
    
    if op == 2: # Válida apenas pro algoritmo específico
        if len(dat.or_v_ex) != len(dat.or_v_in):
            raise TypeError('Vetores originais de extradorso e intradorso devem ter comprimento igual (mesmo grau do polinômio)')
    
    if dat.BPn <= 1:
        raise TypeError('O polinômio de Bernstein deve ter grau de no mínimo 2')
    
    for i in range(dat.cases):
        for j in range(4):
            if dat.coeff_op[i,j] != '!' and dat.coeff_op[i,j] != '^' and dat.coeff_op[i,j] != 'o' and dat.coeff_op[i,j] != 'c' and dat.coeff_op[i,j] != 'k':
                raise TypeError('A Matrix dat.coeff_op deve ser definida com opções !, ^, o, c ou k')
    
    if dat.cases > 1:
        for i in range(1,dat.cases):
            if sum(dat.coeff_op[i,:] == '!') == 4 and dat.coeff_op[0,3] != 'c' and dat.coeff_op[0,3] != 'k':
                raise TypeError('Condição de voo ' + str(i+1) + ' não tem função objetivo definida')
    
    for P in range(dat.cases):
        if dat.coeff_op[P,1] == 'o' and dat.coeff_val[P,1] < 0:
            raise TypeError('Condição de voo ' + str(P+1) + ': CD alvo deve ser maior que zero')
        elif dat.coeff_op[P,1] == 'o' and dat.coeff_val[P,1] == 0:
            dat.coeff_op[P,1] = '^'
            warn('Condição de voo ' + str(P+1) + ': CD = 0 - função objetivo de arrasto trocada de o para ^'),sleep(5) 
    
    if dat.cases > 1 and dat.coeff_op[0,3] == 'c' or dat.cases > 1 and dat.coeff_op[0,3] == 'k':
        if sum(dat.coeff_op[1:,3] != '!') != 0:
            warn('função objetivo de momento c/k: funções objetivo das condições de voo subsequentes serão ignoradas');sleep(5)
            dat.coeff_op[1:dat.cases,3] = repeat('!',dat.cases-1)
    
    return dat