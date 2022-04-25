import numpy as np
import math
# import os
from cosspace_half import cosspace_half
from cosspace import cosspace

def run_cst_TCC2(v_ex,v_in,dat,op=0):

    # Versão nova da função. Agora admite graus diferentes do polinômio.
    # Inputs:
    # - Comprimento da corda c
    # - Grau do polinômio n
    # - Número de pontos NP
    # - Parâmetros da função shape N1 e N2
    # - Vetores com informações do extradorso v_ex
    # - Vetores com informações do intradorso v_in
    # Dos quais c, n, NP, N1 e N2 serão retirados do struct
    # Quando op == 1 as coordenadas serão impressas

    # pegar informações do struct
    c = dat.chord
    n = dat.BPn
    NP = dat.np
    N1 = dat.N1
    N2 = dat.N2
    p_op = dat.p_op

    # Inicializar vetores das coordenadas e pesos
    if p_op == 1:
        x = cosspace(0,c,np)
    else:
        x = cosspace_half(0,c,NP)
    y1 = np.zeros(len(x)); y2 = np.zeros(len(x))
    A1 = np.zeros(n+1); A2 = np.zeros(n+1)

    # Ler informações dos vetores e calcular pesos
    Rle1 = v_ex[0];  Rle2 = v_in[0]
    beta1 = v_ex[n]; beta2 = v_in[n]
    Dz1 = v_ex[n+1]; Dz2 = v_in[n+1]
    A1[0] = np.sqrt(2*Rle1/c)
    A1[1:n] = v_ex[1:n]
    A1[n] = np.tan(np.deg2rad(beta1)) + Dz1/c
    A2[0] = np.sqrt(2*Rle2/c)
    A2[1:n] = v_in[1:n] 
    A2[n] = np.tan(np.deg2rad(beta2)) + Dz2/c

    ## Extradorso ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for p in range(NP): #p = 1:NP
        
        # Calcular o polinômio de Bernstein
        sum = 0
        for r in range(n+1):
            K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
            sum += A1[r]*K*x[p]**r*(1-x[p])**(n-r)
        
        # Calcular a ordenada com as funções class e shape ao mesmo tempo
        y1[p] = x[p]**N1*(1-x[p])**N2*sum + x[p]*Dz1/c

    ## Intradorso ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for p in range(NP):
        
        # Calcular o polinômio de Bernstein
        sum = 0
        for r in range(n+1):
            K = math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
            sum += A2[r]*K*x[p]**r*(1-x[p])**(n-r)
        
        # Calcular a ordenada com as funções class e shape ao mesmo tempo
        y2[p] = -(x[p]**N1*(1-x[p])**N2*sum + x[p]*Dz2/c)

    x = x[...,np.newaxis]
    y1 = y1[...,np.newaxis]
    y2 = y2[...,np.newaxis]
    part1 = np.hstack((np.flip(x,0),np.flip(y1,0)))
    part2 = np.hstack((x[1:NP],y2[1:NP]))
    coo = np.vstack((part1,part2))
    # coo = np.vstack(([np.flip(np.transpose(x)), np.flip(np.transpose(y1))],[np.transpose(x[1:NP]),np.transpose(y2[1:NP])]))
    #figure(1),clf
    #plot(coo(:,1),coo(:,2)),grid on,axis equal,hold on
    #scatter(coo(:,1),coo(:,2))

    # Imprimir coordenadas
    if op == 1:
        ID = open('coordenadas.dat','w')
        # ID = fopen('coordenadas.dat','w');
        with open('coordenadas.dat','w') as f:
            np.savetxt(f, coo, delimiter=' ', newline='\n', header='', footer='', comments='# ')
        
        # fprintf(ID,'#f #f\n',coo');
        ID.close()

    return coo
    
    
