from make_vector_TCC2 import make_vector_TCC2
import numpy as np

def fitness_cst(pop,dat,select2):
    
    for P in range(dat.cases):
    
        # CL (aero[0])
        if dat.coeff_op[P,0] != '!':
            F0 = dat.coeff_F[P,0]
            if dat.coeff_op[P,0] == '^': # Maior CL possível
                temp = make_vector_TCC2(pop,P,0,select2) # Vetor com todos os CLs
                for i in select2:
                    pop[i].score += pop[i].aero[P,0]/max(temp)*F0
            elif dat.coeff_op[P,0] == 'o': # Alcançar um CL específico
                CLtgt = dat.coeff_val[P,0]
                if CLtgt == 0 and dat.aoa[P] == 0 and dat.symm_op == 1: # Se o CL alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
                    for i in select2:
                        pop[i].score += 1 
                elif CLtgt == 0: # Se o CL alvo for igual a zero
                    temp = make_vector_TCC2(pop,P,0,select2) # Vetor com todos os CLs
                    for i in select2:
                        pop[i].score += (1-abs(pop[i].aero[P,0])/max(abs(temp)))*F0
                elif CLtgt > 0: # Se o CL alvo for positivo
                    for i in select2:
                        if pop[i].aero[P,0] >= CLtgt: # Se o CL for maior/igual que o CL alvo
                            pop[i].score += CLtgt/pop[i].aero[P,0]*F0
                        elif pop[i].aero[P,0] < CLtgt: # Se o CL for menor que o CL alvo
                            pop[i].score += pop[i].aero[P,0]/CLtgt*F0
                else: # Se o CL alvo for negativo
                    for i in select2:
                        if pop[i].aero[P,0] <= CLtgt: # Se o CL for menor/igual que o CL alvo
                            pop[i].score += CLtgt/pop[i].aero[P,0]*F0
                        else: # Se o CL for maior que o CL alvo
                            pop[i].score += pop[i].aero[P,0]/CLtgt*F0
        # CD (aero[1])
        if dat.coeff_op[P,1] != '!':
            F1 = dat.coeff_F[P,1]
            if dat.coeff_op[P,1] == '^': # Menor CD possível (tende a zero, mas nunca igual a ou menor que zero)
                temp = make_vector_TCC2(pop,P,1,select2)
                for i in select2:
                    if pop[i].aero[P,1] <= 0:
                        pop[i].score = 0 # Punir aerofólios com resultados irreais
                    else:
                        pop[i].score += (1-pop[i].aero[P,1]/max(temp))*F1
            elif dat.coeff_op[P,1] == 'o': # Alcançar um CD específico
                CDtgt = dat.coeff_val[P,1]
                for i in select2:
                    if pop[i].aero[P,1] >= CDtgt: # Se o CD for maior/igual que o CD alvo
                        pop[i].score += CDtgt/pop[i].aero[P,1]*F1
                    else: # Se o CD for menor que o CD alvo
                        pop[i].score += pop[i].aero[P,1]/CDtgt*F1
        # L/D (aero[2])
        if dat.coeff_op[P,2] != '!':
            F2 = dat.coeff_F[P,2]
            if dat.coeff_op[P,2] == '^': # Maior L/D possível
                temp = make_vector_TCC2(pop,P,2,select2) # Vetor com todos os L/Ds
                for i in select2:
                    pop[i].score += pop[i].aero[P,2]/max(temp)*F2
            elif dat.coeff_op[P,2] == 'o': # Alcançar um L/D específico
                LDtgt = dat.coeff_val[P,2]
                if LDtgt == 0 and dat.aoa[P] == 0 and dat.symm_op == 1: # Se o LD alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
                    for i in select2:
                        pop[i].score += 1
                if LDtgt == 0: # Se o L/D alvo for igual a zero
                    temp = make_vector_TCC2(pop,P,2,select2) # Vetor com todos os L/Ds
                    for i in select2:
                        pop[i].score += (1-abs(pop[i].aero[P,2])/max(abs(temp)))*F2
                elif LDtgt > 0: # Se o L/D alvo for positivo
                    for i in select2:
                        if pop[i].aero[P,2] >= LDtgt: # Se o L/D for maior/igual que o L/D alvo
                            pop[i].score += LDtgt/pop[i].aero[P,2]*F2
                        else:# Se o L/D for menor que o L/D alvo
                            pop[i].score += pop[i].aero[P,2]/LDtgt*F2
                else: # Se o L/D alvo for negativo
                    if pop[i].aero[P,2] <= LDtgt: # Se o L/D for menor/igual que o L/D alvo
                        pop[i].score += LDtgt/pop[i].aero[P,2]*F2
                    else: # Se o L/D for maior que o L/D alvo
                        pop[i].score += pop[i].aero[P,2]/LDtgt*F2 
        # CM (aero[3])
        if dat.coeff_op[P,3] != '!' and dat.coeff_op[P,3] != 'c' and dat.coeff_op[P,3] != 'k':
            F3 = dat.coeff_F[P,3]
            CMtgt = dat.coeff_val[P,3]
            if CMtgt == 0 and dat.aoa[P] == 0 and dat.symm_op == 1: # Se o CM alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
                for i in select2:
                    pop[i].score += 1
            elif CMtgt == 0: # Se o CM alvo for igual a zero
                temp = make_vector_TCC2(pop,P,3,select2) # Vetor com todos os CMs
                for i in select2:
                    pop[i].score += (1-abs(pop[i].aero[P,3])/max(abs(temp)))*F3
            elif CMtgt > 0: # Se o CM alvo for positivo
                for i in select2:
                    if pop[i].aero[P,3] >= CMtgt: # Se o L/D for maior/igual que o L/D alvo
                        pop[i].score += CMtgt/pop[i].aero[P,3]*F3 
                    elif pop[i].aero[P,3] < CMtgt: # Se o L/D for menor que o L/D alvo
                        pop[i].score += pop[i].aero[P,3]/CMtgt*F3 
            else: # Se o CM alvo for negativo
                for i in select2:
                    if pop[i].aero[P,3] <= CMtgt: # Se o CM for menor/igual que o CM alvo
                        pop[i].score += CMtgt/pop[i].aero[P,3]*F3
                    else: # Se o CM for maior que o CM alvo
                        pop[i].score += pop[i].aero[P,3]/CMtgt*F3 
                            
                            
        # Funções objetivas c e k para o coeficiente de momento
        if dat.coeff_op[0,3] == 'c' or dat.coeff_op[0,3] == 'k':
            for i in select2:
                temp = 0
                if dat.coeff_op[0,3] == 'c': # Se o valor for arbitrário
                    CMtgt = np.mean(pop[i].aero[:,3])
                elif dat.coeff_op[0,3] == 'k': # Se o valor for específico
                    CMtgt = dat.coeff_val[0,3]
                
                for P  in range(dat.cases):
                    if CMtgt == 0: # Se o CM alvo for igual a zero
                        temp2 = make_vector_TCC2(pop,P,3,select2) # Vetor com todos os CMs
                        temp += (1-abs(pop[i].aero[P,3])/np.max(abs(temp2)))
                    elif CMtgt > 0: # Se o CM alvo for positivo
                        if pop[i].aero[P,3] >= CMtgt: # Se o CM for maior/igual que o CM alvo
                            temp += CMtgt/pop[i].aero[P,3]
                        else: # Se o CM for menor que o CM alvo
                            temp += pop[i].aero[P,3]/CMtgt
                        
                    else: # Se o CM alvo for negativo
                        if pop[i].aero[P,3] <= CMtgt: # Se o CM for menor/igual que o CM alvo
                            temp += CMtgt/pop[i].aero[P,3]
                        else: # Se o CM for maior que o CM alvo
                            temp += pop[i].aero[P,3]/CMtgt
                
                pop[i].score += (temp/dat.cases)*dat.coeff_F[0,3]
                
                
    return pop