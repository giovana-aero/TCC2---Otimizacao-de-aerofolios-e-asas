import numpy as np

def cosspace_half(startP,endP,N=100):

    # Versão modificada da funçãoo cosspace. Gera apenas metade da distribuição.
    x = np.zeros(N); x[N-1] = endP
    angleInc = np.pi/(N-1)/2

    curAngle = angleInc
    for i in range(1,N-1):
        x[i]=endP*(1-np.cos(curAngle))
        curAngle += angleInc

    return x