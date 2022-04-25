clc,clear


% Objetivo: testar a função de geração de malha pro apame

% Dados da planta
pop.type = 1; % 0 -> Trapezoidal simples, 1 -> Trapezoidal dupla
pop.b = 20;
pop.b1 = 10;
pop.c_r = 5;
pop.c_m = 2;
pop.c_t = 1;
pop.tw_m = 'L';
pop.tw_t = 0;

% Partes do struct que serão preenchidas posteriormente
pop.coo_r = [];
pop.coo_m = [];
pop.coo_t = [];
pop.NODE = [];
pop.ELEM = [];
pop.sec_N = [];
pop.aero = [];

% Aerofólios CST
pop.v_ex_r = [0.0100, 0.1300, 0.2900, 0.2300, 16.0000, 0.0000];
pop.v_in_r = [0.0100, 0.0800, 0.0800, 0.0500, 13.0000, 0.0000];
pop.v_ex_m = [0.0100, 0.1600, 0.3000, 0.1100, 21.0000, 0.0000];
pop.v_in_m = [0.0100, 0.0200, 0.0200, 0.0200, -1.0000, 0.0000];
pop.v_ex_t = [0.0100, 0.2500, 0.1800, 0.2900, 19.0000, 0.0000];
pop.v_in_t = [0.0100, 0.0700, 0.0500, 0.0300, 1.0000, 0.0000];
dat.BPn = length(pop.v_ex_r)-2;
dat.N1 = 0.5; dat.N2 = 1;

% Número de pontos
dat.np = 30;
dat.p_op = 1;
dat.chord = 1;

% Dados da malha de painéis
dat.nb = [1,1]; % Número de seções intermediárias (raiz/ponta) [número de seções,desligado (0)] ou [seções por metro,ligado(1)]
dat.nb1 = 'L'; % Número de seções intermediárias (raiz/meio)
dat.nb2 = []; % Número de seções intermediárias (meio/ponta)

% Dados da simulação
dat.cases = 1;
dat.v_ref = 100; % [m/s]
dat.rho = 1.225; % [kg/m^3]
dat.p_atm = 101325; % [Pa]
dat.aoa = 0; % [graus]
dat.mach = 0.; 

% Obter coordenadas dos aerofólios
%pop.af_r = [str2num(pop.naca_r(1)),str2num(pop.naca_r(2)),str2num(pop.naca_r(3:4))]; % Raiz
%if pop.type == 1
%    pop.af_m = [str2num(pop.naca_m(1)),str2num(pop.naca_m(2)),str2num(pop.naca_m(3:4))]; % Meio
%end    
%pop.af_t = [str2num(pop.naca_t(1)),str2num(pop.naca_t(2)),str2num(pop.naca_t(3:4))]; % Ponta

% Obter geometria (malha) do APAME
figure(1),clf,figure(2),clf
pop = run_apame_mesher_cst(pop,dat,2,2,1,2,'Título');

% Escrever o arquivo de entrada e fazer a simulação
pop.aero = run_apame(pop,dat);

% Imprimir valores
if pop.aero ~= 'n'
    fprintf('CL  = %.8f\nCD  = %.8f\nL/D = %.8f\nCM  = %.8f\n',pop.aero(1),pop.aero(2),pop.aero(3),pop.aero(4))
else
    disp('A simulação não convergiu')
end

