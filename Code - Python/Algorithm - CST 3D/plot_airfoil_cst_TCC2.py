import matplotlib.pyplot as plt
from numpy import genfromtxt
import locale
locale.setlocale(locale.LC_ALL, 'pt_BR.UTF-8') # Trocar separadores decimais
plt.ticklabel_format( axis='both', style='', scilimits=None, useOffset=None, useLocale=True, useMathText=None)

# - As duas primeiras entradas (coordenadas,op) são obrigatórias para todos os casos
# - A terceira entrada é obrigatória para op = 1 e op = 2
# - A quarta entrada é obrigatória para op = 2

def plot_airfoil_cst_TCC2(op=0,coordenadas=None,loop=None,pop=None):
    
    if op == 0: # Ler coordenadas a partir de um arquivo de texto
        coo = genfromtxt('coordenadas.dat',delimiter=' ')
        plt.plot(coo[:,0],coo[:,1]) #CONSERTAR ISTO
    
    elif op == 1: # Pôr no título apenas o número da iteração
        plt.plot(coordenadas[:,0],coordenadas[:,1])
        plt.title('Iteração ' + str(loop))
        
    elif op == 2: # Pôr no título as informações aerodinâmicas
        plt.plot(coordenadas[:,0],coordenadas[:,1])
        aero = pop.aero[0,:]
        plt.title('Iteração ' + str(loop) + ': CL = ' + str(aero[0]) + ', CD = ' + str(aero[1]) + ', L/D = ' + str(aero[2]) + ', CM = ' + str(aero[3]))
        
    elif op == 3: # Fazer com estilo tracejado para fins de comparação
        plt.plot(coordenadas[:,0],coordenadas[:,1],'--')
        
    elif op == 4: # Fazer um gráfico sem nada no título
        plt.plot(coordenadas[:,0],coordenadas[:,1])