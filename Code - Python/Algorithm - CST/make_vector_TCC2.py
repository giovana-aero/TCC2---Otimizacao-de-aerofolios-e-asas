import numpy as np

def make_vector_TCC2(pop,P,op,select2):
    # Esta função toma as informações das matrizes aero e transforma em vetores
    # Atualização: especificação da linha desejada da matriz aero
    
    v = np.zeros(len(pop))
    for i in select2:
        v[i] = pop[i].aero[P,op]
    
    return v