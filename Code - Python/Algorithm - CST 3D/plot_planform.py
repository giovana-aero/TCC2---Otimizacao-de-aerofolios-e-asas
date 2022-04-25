import math
import numpy as np
import matplotlib.pyplot as plt

def plot_planform(pop,color):
    
    c_r = pop.c_r
    c_t = pop.c_t
    tw_t = pop.tw_t
    b = pop.b
    
    if pop.type == 0: # Se a asa for trapezoidal simples
        
        x = [0,(c_r-c_t*np.cos(math.radians(tw_t)))/2,(c_r + c_t*np.cos(math.radians(tw_t)))/2,c_r]
        x = np.hstack((x,np.flip(x[0:-1])))
        y = [0,b/2,b/2,0]
        y = np.hstack((y,-1*np.flip(y[0:-1])))
        
        plt.plot(x,y,color)
        
    elif pop.type == 1: # Se a asa for trapezoidal dupla
        
        b1 = pop.b1
        c_m = pop.c_m
        tw_m = pop.tw_m
        # Considerar opção 'L' da corda do meio
        if c_m == 'L':
            c_m = 2*(c_t-c_r)/b*b1/2 + c_r
        
        # Considerar opção 'L' da torção do meio
        if tw_m == 'L' and not isinstance(tw_t,str): # Definir a torção do meio em termos da torção da ponta
            tw_m = tw_t*b1/b
        elif tw_t == 'L' and not isinstance(tw_m,str): # Definir a torção da ponta em termos da torção do meio
            tw_t = tw_m*b/b1
        
        x = [0,(c_r - c_m*np.cos(math.radians(tw_m)))/2,(c_r - c_t*np.cos(math.radians(tw_t)))/2,(c_r + c_t*np.cos(math.radians(tw_t)))/2,(c_r + c_m*np.cos(math.radians(tw_m)))/2,c_r]
        x = np.hstack((x,np.flip(x[0:-1])))
        y = [0,b1/2,b/2,b/2,b1/2,0]
        y = np.hstack((y,-1*np.flip(y[0:-1])))
        plt.plot(x,y,color)