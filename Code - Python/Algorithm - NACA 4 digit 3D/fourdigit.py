import numpy as np
import math
from cosspace import cosspace
from cosspace_half import cosspace_half

def fourdigit(naca_num,n_p=100,x_op=1):
    
    # Entradas:
    # Nome do perfil (vetor 1x3: [m,p,t])
    # np - número de pontos (abscissa)
    # x_op - opção de geração de pontos (1 pra cosspace, outro valor pra cosspace_half)
    
    #if nargin == 4,x_op=1;end
    
    if x_op == 1:  # Opção padrão
        x = cosspace(0,1,n_p)
    else:
        x = cosspace_half(0,1,n_p)
    
    
    m = naca_num[0]/100
    p = naca_num[1]/10
    t = naca_num[2]/100
    
    # Distribuição de espessura
    # (o último valor muda de -0.1015 para -0.1036 no caso de bordo de fuga
    # fechado)
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036
    yt = 5*t*(a0*x**0.5+a1*x+a2*x**2+a3*x**3+a4*x**4)
    
    if m == 0 and p == 0: # Perfil simétrico
        xU = x
        yU = yt
        xL = x
        yL = -yt
        
    else:  # Perfil com curvatura
        
        # Distribuição de curvatura
        yc = np.zeros(len(x))
        for i in range(len(x)):
            if x[i] >= 0 and x[i] <= p:
                yc[i] = m/p**2*(2*p*x[i]-x[i]**2)
            elif x[i] >= p and x[i] <= 1:
                yc[i] = (m/(1-p)**2)*((1-2*p)+2*p*x[i]-x[i]**2)
        
        derivada = np.zeros(len(x))
        for i in range(len(x)):
            if x[i] >= 0 and x[i] <= p:
                derivada[i] = 2*m/p**2*(p-x[i])
            elif x[i] >= p and x[i] <= 1:
                derivada[i] = 2*m/(1-p)**2*(p-x[i])
        
        theta = math.atan(derivada)
        xU = x - yt*math.sin(theta)
        yU = yc + yt*math.cos(theta)
        xL = x + yt*math.sin(theta)
        yL = yc - yt*math.cos(theta)
        
    
    # Coordenadas
    coo = np.vstack((np.hstack((np.flip(xU),np.flip(yU))),np.hstack((xL[1:],yL[1:]))))
    # coo2 = 
    # coo = [flip(xU'),flip(yU');
    # xL(2:end)',yL(2:end)'];
           
    return coo