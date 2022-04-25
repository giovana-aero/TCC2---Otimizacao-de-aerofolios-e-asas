from ypstruct import structure
from run_apame_mesher_naca4 import run_apame_mesher_naca4
from run_apame import run_apame
import matplotlib.pyplot as plt


# Dados da asa
pop = structure()
pop.type = 0;
pop.b = 20; # Envergadura
pop.b1 = 14;
pop.c_r = 5; # Corda da raiz
pop.c_m = 4; 
pop.c_t = 1; # Corda da ponta (DEVE ser menor que a da raiz)
pop.tw_m = 'L';
pop.tw_t = -0; # Torção geométrica na ponta
pop.af_r = [2,4,12];
pop.af_m = [2,4,12];
pop.af_t = [0,0,10];

# (enflechamentos são sempre referentes ao bordo de ataque da asa)
# Opção 'Z' aqui calcula o enflechamento automaticamente de modo que o enflechamento da linha c/2 seja sempre zero
pop.sweep = 'Z'; # enflechamento da asa (trapezoidal simples) [graus]
pop.sweep1 = 'Z'; # enflechamento da primeira seção (trapezoidal dupla) [graus]
pop.sweep2 = 'Z'; # enflechamento da segunda seção (trapezoidal dupla) [graus]


dat = structure()
dat.np = 30; # Número de pontos na geração de ordenadas nos aerofólios
dat.np_op = 1; # 1 -> cosspace, 0 -> cosspace_half
dat.nb = [1,1];
dat.nb1 = 'L';
dat.nb2 = 2;
# Dados da simulação
dat.cases = 1
dat.v_ref = [100,100,100]; # Velocidades de referência [m/s] 
dat.rho = [1.225,1.225,1.225]; # Densidades do ar [kg/m^3] 
dat.p_atm = [101325,101325,101325]; # Pressões do ar [Pa] (irrelevante neste algoritmo)
dat.mach = [0,0.,0.2]; # Números de Mach
dat.reynolds = [1e6,1e6,1e6];           # Valores dos números de Reynolds para as simulações (irrelevante neste algoritmo))
dat.aoa = [0,2,4];    



# figure(1),clf#,figure(2),clf
pop = run_apame_mesher_naca4(pop,dat,1)
#pop = run_apame_mesher_naca4_uuuh(pop,dat,2,0,1,2,'teste');

# a,b = run_apame(pop,dat)
aero = run_apame(pop,dat)
print(aero)

# fid = open('apame_input.inp','w')
# fid.write('# GEOMETRY\n\n')
# fid.write('NODES #d\n',pop.NODE.shape[0])
# fid.write('#.10f #.10f #.10f\n',pop.NODE)
# # fid.write('\nPANELS #d\n',pop.PANEL.shape[0])
# # fid.write('#d #d #d #d #d #d #d #d #d\n',pop.PANEL[1:-pop.sec_N+1,:])
# # fid.write('#d #d #d #d #d #d #d\n',pop.PANEL[-pop.sec_N+2:,0:6])
# fid.close()

