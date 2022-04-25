import numpy as np

def make_vector_TCC2_3D(pop,P,coef,op,select2,v_ref=None,rho=None):
    
    # Esta função transforma informações dos vetores aero e transforma em
    # um vetor
    # Atualização: especificação da linha desejada da matriz aero
    
    # op = 1 -> coeficientes 
    # op = 2 -> forças de sustentação e arrasto
    # op = 3 -> momento de arfagem
    
    v = np.zeros(len(pop))
    if op == 1:
        for i in select2:
            v[i] = pop[i].aero[P,coef]
        
    elif op == 2:
        for i in select2:
            v[i] = pop[i].aero[P,coef]*1/2*rho[P]*v_ref[P]**2*pop[i].S;
        
    else:
        for i in select2:
            v[i] = pop[i].aero[P,coef]*1/2*rho[P]*v_ref[P]**2*pop[i].S*pop[i].mac
    
    return v