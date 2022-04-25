import numpy as np
import matplotlib.pyplot as plt
import locale
locale.setlocale(locale.LC_ALL, 'pt_BR.UTF-8') # Trocar separadores decimais
plt.ticklabel_format( axis='both', style='', scilimits=None, useOffset=None, useLocale=True, useMathText=None)

def plot_airfoil_naca4_TCC2(pop,loop):
    
    
    Chord = 1
    x = np.arange(0,Chord,0.01)
    m = pop.m/100
    p = pop.p*Chord/10
    t = pop.t/100
    #aero = pop.aero

    if m == 0:
        Symm = 1
    else:
        Symm = 0

    if Symm == 1:
        y_upper = 5*t*Chord*(0.2969*np.sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord)**2+0.2843*(x/Chord)**3-0.1015*(x/Chord)**4)
        y_lower = -y_upper
        x_upper = x
        x_lower = x
    else:
        y_camber = np.zeros(len(x))
        dy_camber = np.zeros(len(x))
        for i in range (len(x)):
            if x[i]/Chord <= p:
                y_camber[i] = m*x[i]/p**2*(2*p-x[i]/Chord)
                dy_camber[i] = 2*m/p**2*(p-x[i]/Chord)
            else:
                y_camber[i] = m*(Chord-x[i])/(1-p)**2*(1+x[i]/Chord-2*p)
                dy_camber[i] = 2*m/(1-p)**2*(p-x[i]/Chord)

        y_t = 5*t*Chord*(0.2969*np.sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord)**2+0.2843*(x/Chord)**3-0.1015*(x/Chord)**4);
        theta = np.arctan(dy_camber)
        x_upper = x-y_t*np.sin(theta)
        x_lower = x+y_t*np.sin(theta)
        y_upper = y_camber+y_t*np.cos(theta)
        y_lower = y_camber-y_t*np.cos(theta)

    # Traçar o perfil
    plt.plot(np.concatenate((np.flip(x_upper),x_lower)),np.concatenate((np.flip(y_upper),y_lower)))
    #plot([flip(x_upper) x_lower],[flip(y_upper) y_lower])
    #axis equal,grid on
    plt.grid(True);plt.axis('equal')
    plt.title('Iteração ' + str(loop) + ': NACA ' + str(pop.m) + str(pop.p) + str(pop.t))
    #title(['Iteração ' num2str(loop) ': NACA ' num2str(pop.m) num2str(pop.p) num2str(pop.t)])

    #xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
    #new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
    #set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)


    
