from math import radians,cos,sin
from numpy import array,matmul

def airfoil_rotation(coo,th):

    # Centro de rotação (sempre no ponto de quarto de corda como padrão)
    c = 1/4
    
    # Fazer a rotação
    th = -th # As convenções de sinal (matemática e aeronáutica) são inversas, portanto, realiza-se a troca aqui
    
    
    R = array(([cos(radians(th)),-sin(radians(th))],
                   [sin(radians(th)),cos(radians(th))]))
    coo[:,0] -= c  # Transladar as coordenadas pra que a rotação ocorra em torno do ponto desejado
    # coo_R = (R*coo')';
    coo_R = matmul(R,coo.transpose()).transpose()
    
    # (R*coo.transpose()).transpose
    coo_R[:,0] += c     # Mover a geometria de volta à posição original
    
    
    
    return coo_R