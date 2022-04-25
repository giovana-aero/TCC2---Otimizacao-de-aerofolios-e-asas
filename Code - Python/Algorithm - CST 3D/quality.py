import numpy as np

def quality(coo,dat):
# Checagem de qualidade dos perfis CST

    NP = dat.np
    c = dat.chord

    # Separar as superfícies
    EX = np.flip(coo[0:NP,1],0) 
    IN = coo[(NP-1):(NP*2-1),1]
    # Definir a espessura
    thickness = EX - IN

    check = 1
    while True:
        
        # 1 - Checar se há interseção
        if sum(thickness[1:-1]<=0) > 0:
            check = 0;break
        
        # 2 - Checar as inclinações do extradorso
        # A inclinação do extradorso não pode mudar duas ou mais vezes
        slope = np.zeros(len(EX))
        for i in range(len(EX)-1):
            slope[i] = EX[i+1] - EX[i]
        v = slope > 0
        
        # Isto detecta quantas mudanças de inclinações existem
        counter = 0; num = 0
        for i in range (len(v)):
            if v[i] == num:
                counter += 1
                if counter % 2 == 0:
                    num = 0
                else:
                    num = 1

        if counter >= 2:
            check = 0; break
        
        # 3 - Checar as inclinações do intradorso
        # A inclinação do intradorso não pode mudar três ou mais vezes
        slope = np.zeros(len(IN))
        for i in range(len(IN)-1):
            slope[i] = IN[i+1] - IN[i]
        v = slope > 0
        
        counter = 0; num = 0
        for i in range(len(v)-1):
            if v[i] == num:
                counter += 1
                if counter % 2 == 0:
                    num = 0
                else:
                    num = 1
        
        if counter >= 3:
            check = 0; break
        
        # 4 - Garantir que o perfil não seja fino demais
        if np.mean(thickness)/c < 0.04:
            check = 0; break
        
        break

    return check