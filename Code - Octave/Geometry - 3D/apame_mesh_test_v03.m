% Objetivo: montar a malha de uma asa pro APAME
clc,clear,fclose('all');

% Características 
% - Asa trapezoidal dupla
% - Um perfil pra raiz, um pro meio e outro pra ponta
% - Pode-se especificar o comprimento das cordas
% - Pode-se aplicar torção à ponta da asa
% Configurações de torção geométrica
% - torções independentes: torção da ponta e do meio são valores determinados pelo usuário;
% - torções lineares: a torção do meio é calculada a partir da torção da ponta ou vice-versa

% Nota sobre as coordenadas de aerofólios: este script trabalha de modo que seja
% necessário que os aerofólios carregados tenham sempre o mesmo número de pontos.
% Isso não é um problema no contexto do algoritmo de otimização. No entanto, se o
% o intuito for gerar uma malha em outra aplicação, é necessário garantir esse
% requisito de números iguais de pontos. Para isso, pode-se interpolar as coordenadas
% em mãos, ou inserí-las no CST reverso e gerar o aerofólio com a quantidade desejada
% de pontos




% NOTA: essa funcionalidade da especificação das torções deve ser aplicada ao código
% do python com a API do OpenVSP


% Dados da asa
b = 10; % Envergadura (total)
b1 = 5; % Envergadura da primeira metade da asa (raiz ao meio)
c_r = 5; % Corda da raiz
c_m = 2; % Corda do meio 
c_t = 1; % Corda da ponta (DEVE ser menor que a da raiz)
tw_m = 'L'; % Torção geométrica no meio 
tw_t = -5; % Torção geométrica na ponta
% Em qualquer uma das torções, pode-se especificar a opção 'L'
% Por exemplo: se tw_m = 'L', a torção é aplicada à asa completa, linearmente, 
% em termos da torção da ponta
% NOTA: tw_t <= tw_m
% Para a corda do meio: a opção 'L' essencialmente torna a asa em uma trapezoidal
% simples (a graça é ter mais opções para alterar o aerofólio)


% Aerofólios CST
v_ex_r = [0.04,.2,.3,.2,20,0]; v_in_r = [.01,.1,.1,.2,0,0]; 
v_ex_m = [0.04,.2,.3,.1,20,0]; v_in_m = [.01,.1,.1,.2,-10,0];
v_ex_t = [0.01,.1,.1,.1,10,0]; v_in_t = [.01,.1,.1,.2,10,0]; 
dat.chord = 1;
dat.BPn = 4;
dat.np = 11;
dat.N1 = 0.5;
dat.N2 = 1;
dat.p_op = 0;
% Carregar coordenadas dos aerofólios. Como o contorno é fechado, ignora-se o último par de coordenadas
coo_r = run_cst_TCC2(v_ex_r,v_in_r,dat); coo_r = coo_r(1:end-1,:)*c_r;
if c_m == 'L',c_m = c_t*b1/b;end
coo_m = run_cst_TCC2(v_ex_m,v_in_m,dat); coo_m = coo_m(1:end-1,:)*c_m;
coo_t = run_cst_TCC2(v_ex_t,v_in_t,dat); coo_t = coo_t(1:end-1,:)*c_t;

% Nota: aerofólios *devem* ter o bordo de fuga fechado

% Dados da malha de painéis
far = b*2; % Comprimento dos painéis de trilha (a partir do bordo de fuga)
b2 = b - b1;
nb1 = 5;    % Número de seções intermediárias entre meio/raiz (considerando apenas um lado da asa)
nb2 = 6;    % Número de seções intermediárias entre ponta/meio 
% NOTA: por hora, o código só funciona como devido quando nb1=nb2 (corrigir isto depois)

% Configurações da simulação
sim_op = 1; % Fazer a simulação no apame ao final?
v_ref = 100; % [m/s]
aoa = 0; % [graus]
rho = 1.225; % [kg/m^3]
p_atm = 101325; % [Pa]

% Alguns dados a mais
sec_af_N = size(coo_r,1); % Número de nós por seção (cada aerofólio)
sec_N = 5 + nb1*2 + nb2*2; % Número de seções transversais
    
% Distribuição de torções geométricas
if tw_m == 'L' && ~ischar(tw_t) % Definir a torção do meio em termos da torção da ponta
    tw_m = tw_t*b1/b;
    tw_V = linspace(0,1,3+nb1+nb2)*tw_t;
    tw_V1 = flip(tw_V(2:nb1+1));
    tw_V2 = flip(tw_V(3+nb1:end-1));
elseif tw_t == 'L' && ~ischar(tw_m) % Definir a torção da ponta em termos da torção do meio
    tw_t = tw_m*b/b1;
    tw_V = linspace(0,1,3+nb1+nb2)*tw_t;
    tw_V1 = flip(tw_V(2:nb1+1));
    tw_V2 = flip(tw_V(3+nb1:end-1));
else % Aplicar ambos os valores de torção especificados
    tw_V1 = linspace(0,1,2+nb1)*tw_m;
    tw_V2 = linspace(tw_m/tw_t,1,2+nb2)*tw_t;
    tw_V1 = flip(tw_V1(2:end-1));
    tw_V2 = flip(tw_V2(2:end-1));
end
% (Em todos os casos acima, as torções intermediárias são definidas linearmente)

% Rotacionar o perfil da ponta
coo_t_R = airfoil_rotation(coo_t,tw_t);

% Obter nós ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NODE = zeros(sec_af_N*(1+nb2),3);
dx_t = (c_r - c_t*cosd(tw_t))/2; % Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
%tw_m = tw_t*(b1/2)/(b/2);
dx_m = (c_r - c_m*cosd(tw_m))/2;
NODE(1:sec_af_N,:) = [coo_t_R(:,1)+dx_t,repmat(-b/2,sec_af_N,1),coo_t_R(:,2)]; % Ponta esquerda da asa
%NODE = [coo_t(:,1),repmat(-b/2,sec_af_N,1),coo_t(:,2)]; % Ponta esquerda da asa

%figure(3),clf
%plot(coo_t(:,1),coo_t(:,2))
%axis equal,hold on,grid on
%plot(coo_t_R(:,1),coo_t_R(:,2))
%plot(NODE(1:sec_af_N,1),NODE(1:sec_af_N,3))

% Gerar seções intermediárias (ponta/meio)
if nb2 > 0
    % Criar struct que guarda seções de asa intermediárias
    % (isto será um auxílio devido à natureza simétrica a asa)
    wing_sec2.coo = []; wing_sec2.dx = 0;
    wing_sec2 = repmat(wing_sec2,nb2,1); % Inicializar
    op_vec = linspace(1,0,2+nb2); op_vec = op_vec(2:end-1); % Definição do formato da interpolação em função dos originais
    dx_t2 = (c_m - c_t*cosd(tw_t))/2; % Mesma coisa que dx_t, mas em relação à seção do meio
    
    % Encontrar as coordenadas das seções (fazer interpolações)
    for i = 1:length(op_vec)
        % Fazer interpolação e aplicar torção geométrica
        tw_sec = tw_V2(i);%*cosd(c_t*op_vec(i)+c_r*(1-op_vec(i)));
        wing_sec2(i).coo = airfoil_interpolation(coo_m,coo_t,op_vec(i),-op_vec(i)*((b-b1)/2)-b1/2,tw_sec);
        % Aplicar o translado
        wing_sec2(i).dx = op_vec(i)*dx_t2 + dx_m; 
        wing_sec2(i).coo(:,1) = wing_sec2(i).coo(:,1) + wing_sec2(i).dx;
        % Adicionar ao lado esquerdo da asa
        NODE(sec_af_N*i+1:sec_af_N*(i+1),:) = wing_sec2(i).coo;
    end
        
end

% Rotacionar perfil do meio
coo_m_R = airfoil_rotation(coo_m,tw_m);

% Adicionar coordenadas do perfil do meio
NODE = [NODE;coo_m_R(:,1)+dx_m,zeros(sec_af_N,1)-b1/2,coo_m_R(:,2)];
%NODE(1:sec_af_N,:) = [coo_t_R(:,1)+dx_t,repmat(-b/2,sec_af_N,1),coo_t_R(:,2)]

% Gerar seções intermediárias (meio/raiz)
if nb1 > 0
    temp = zeros(sec_af_N*(nb1),3);
    % Criar struct que guarda seções de asa intermediárias
    % (isto será um auxílio devido à natureza simétrica a asa)
    wing_sec1.coo = []; wing_sec1.dx = 0;
    wing_sec1 = repmat(wing_sec1,nb1,1); % Inicializar
    op_vec = linspace(1,0,2+nb1); op_vec = op_vec(2:end-1); % Definição do formato da interpolação em função dos originais
    
    % Encontrar as coordenadas das seções (fazer interpolações)
    for i = 1:length(op_vec)
        % Fazer interpolação e aplicar torção geométrica
        tw_sec = tw_V1(i);%*cosd(c_t*op_vec(i)+c_r*(1-op_vec(i)));
        wing_sec1(i).coo = airfoil_interpolation(coo_r,coo_m,op_vec(i),-op_vec(i)*(b1/2),tw_sec);
        % Aplicar o translado
        wing_sec1(i).dx = op_vec(i)*dx_m; % op_vec(i)*(b/2)*dx_t/(b/2)
        wing_sec1(i).coo(:,1) = wing_sec1(i).coo(:,1) + wing_sec1(i).dx;
        % Adicionar ao lado esquerdo da asa
        temp(sec_af_N*(i-1)+1:sec_af_N*i,:) = wing_sec1(i).coo;
    end
    NODE = [NODE;temp];    
end

% Adicionar coordenadas do perfil da raiz
NODE = [NODE;coo_r(:,1),zeros(sec_af_N,1),coo_r(:,2)];

% Adicionar seções intermediárias ao lado direito da asa (raiz/meio)
if nb1 > 0
%    NODE = [NODE;zeros(sec_af_N*nb1,3)];
    temp = zeros(sec_af_N*nb1,3);
    k = 1;
    for i = length(op_vec):-1:1
        wing_sec1(k).coo(:,2) = -wing_sec1(k).coo(:,2); % Inverter o sinal da coordenada y das seções intermediárias
        temp(sec_af_N*i-sec_af_N+1:sec_af_N*i,:) = wing_sec1(k).coo;
        k = k + 1;
    end
    NODE = [NODE;temp];
end

% Adicionar perfil do meio do lado direito da asa
NODE = [NODE;coo_m_R(:,1)+dx_m,zeros(sec_af_N,1)+b1/2,coo_m_R(:,2)];

% Adicionar seções intermediárias ao lado direito da asa (meio/ponta)
if nb2 > 0
%    NODE = [NODE;zeros(sec_af_N*nb,3)];
    temp = zeros(sec_af_N*nb2,3);
    k = 1;
    for i = length(op_vec):-1:1
        wing_sec2(k).coo(:,2) = -wing_sec2(k).coo(:,2); % Inverter o sinal da coordenada y das seções intermediárias
        temp(sec_af_N*i-sec_af_N+1:sec_af_N*i,:) = wing_sec2(k).coo;
        k = k + 1;
    end
    NODE = [NODE;temp];
end

% Adicionar ponta direita
NODE = [NODE;coo_t_R(:,1)+dx_t,zeros(sec_af_N,1)+b/2,coo_t_R(:,2)];

% Teste: fazer um gráfico das seções da asa (raiz à ponta)
figure(2),clf,grid on,axis equal,hold on
for i = 1:(sec_N-1)
%    figure(2),clf,grid on,axis equal
    plot(NODE(sec_af_N*i-sec_af_N+1:sec_af_N*i,1),NODE(sec_af_N*i-sec_af_N+1:sec_af_N*i,3))
%    input('')
end

% Fazer um gráfico da planta da asa
figure(3),clf
scatter(NODE(:,1),NODE(:,2)),grid on,axis equal




% Teste: fazer um gráfico dos perfis carregados
%figure(1),clf
%plot(coo_r(:,1),coo_r(:,2)),axis equal,grid on,hold on
%plot(coo_t(:,1),coo_t(:,2))


% Gerar painéis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Número de painéis: (sec_af_N*(3 + nb*2 - 1)) + ((3 + nb*2 - 1))
%               (painéis na superfície das asas) + (painéis de trilha)
panel_surf = zeros(sec_af_N*(sec_N-1),9); % Inicializar painéis da superfície
% Devido à natureza das geometrias, não haverão elementos triangulares 
% Nota: painéis da superfície que se localizam no bordo de fuga não tem um dos
% painéis adjacentes. Nesse caso, insere-se o valor 0

% Montar os elementos da superfície
v = 1:(sec_N-1);
k = 1;
for i = v
    for j = 1:sec_af_N-1
        panel_surf(k,1:5) = [1,k,sec_af_N+k,sec_af_N+k+1,k+1];
%        j+sec_af_N*(i-1)
        k = k + 1;
    end    
    panel_surf(k,1:5) = [1,k,sec_af_N+k,sec_af_N*i+1,sec_af_N*(i-1)+1];
    k = k + 1;
end
% Adicionar numeração de painéis adjacentes
panel_surf(1,6:end) = [sec_af_N+1,2,0,0]; % Primeiro painel
panel_surf(end,6:end) = [sec_af_N*(sec_N-2),sec_af_N*(sec_N-1)-1,0,0]; % Último painel
for i = 2:size(panel_surf,1)-1
    if i <= sec_af_N % Ponta esquerda da asa
        panel_surf(i,6:end) = [i-1,i+sec_af_N,i+1,0];
    elseif i > size(panel_surf,1)-sec_af_N+1 % Ponta direita da asa
        panel_surf(i,6:end) = [i-1,0,i+1,i-sec_af_N];
    elseif i == sec_af_N % Painel do intradorso da ponta esquerda (bordo de fuga)
        panel_surf(i,6:end) = [sec_af_N*2,0,0,sec_n-1];
    elseif i == sec_af_N*(1+2*(nb1+nb2))+1 % Painel do extradorso da ponta direita(bordo de fuga)
        panel_surf(i,6:end) = [0,0,sec_af_N*(1+(nb1+nb2)*2)+2,sec_af_N*(1+(nb1+nb2))+1];
    else % Todos os outros pontos
        panel_surf(i,6:end) = [i-1,i+sec_af_N,i+1,i-sec_af_N];
    end
end

% Adicionar os nós da trilha da asa
% Cada um será posicionado diretamente atrás de sua respectiva seção de asa a uma distância far
far_nodes = zeros(sec_N-1,3);
for i = 1:sec_N
    far_nodes(i,:) = [NODE(sec_af_N*i-sec_af_N+1,1)+far,NODE(sec_af_N*i-sec_af_N+1,2:3)];
end
S = size(NODE,1);
NODE = [NODE;far_nodes];

% Teste: fazer um gráfico dos nós
figure(1),clf
scatter3(NODE(:,1),NODE(:,2),NODE(:,3)),axis equal
xlabel('x'),ylabel('y'),zlabel('z')
% Enumerar os nós
%for i = 1:size(NODE,1)   
%%    text(x(n),y(n),num2str(n))
%    text(NODE(i,1),NODE(i,2),NODE(i,3),num2str(i))
%end

% Montar os painéis da trilha
%num_sec = 3 + 2*nb;
panel_far = zeros(sec_N-1,9);
%S = size(panel_surf,1);
k = 1;
for i = 1:sec_N-1
    panel_far(i,:) = [10,k,S+i,S+i+1,k+sec_af_N,k,k+sec_af_N-1,0,0];
    k = k + sec_af_N;
end
PANEL = [panel_surf;panel_far];

%for i = 1:(3 + nb*2)
%    num = sec_af_N*i-sec_af_N+1;
%    disp(NODE(num,:))
%    disp(num)
%end

%op_vec = linspace(1,0,2+nb);
%clc,for i = 1:(nb+2)
%    disp(sec_af_N*i-sec_af_N+1)
%end

% Imprimir o arquivo de entrada ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

fid = fopen('wings_rewrite_v03.inp','w');
fprintf(fid,'APAME input file\nVERSION 3.1\n\n\n');

fprintf(fid,'# FLOW PARAMETERS\n\n');
fprintf(fid,'# Airspeed [m/s]\n');
fprintf(fid,'AIRSPEED %f\n',v_ref);
fprintf(fid,'# Air density\n');
fprintf(fid,'DENSITY %f\n',rho);
fprintf(fid,'# Atmospheric pressure\n');
fprintf(fid,'PRESSURE %f\n',p_atm);
fprintf(fid,'# Prandtl-glauert corretion:\n#  0-no correction\n#  *-Mach number\n');
fprintf(fid,'MACH 0\n');
fprintf(fid,'# Number of cases\n# Angle of attack [degrees]\n# sideslip angle [degrees]\n');
fprintf(fid,'CASE_NUM 1\n')
fprintf(fid,'%f\n0\n\n\n',aoa);

fprintf(fid,'# REFERENCE VALUES\n\n');
fprintf(fid,'# Wingspan [m]\n');
fprintf(fid,'WINGSPAN %f\n',b);
fprintf(fid,'# Mean aerodynamic chord [m]\n');
fprintf(fid,'MAC %f\n',(c_r + c_m + c_t)/3);
fprintf(fid,'# Wing surface area [m^2]\n');
fprintf(fid,'SURFACE %f\n',b);
fprintf(fid,'# Reference point [m]\n');
fprintf(fid,'ORIGIN *\n0 0 0\n\n\n');

fprintf(fid,'# SOLVER PARAMETERS\n\n');
fprintf(fid,'# Singularity method:\n#  0-constant source/doublet\n#  1-constant doublet\n');
fprintf(fid,'METHOD 0\n');
fprintf(fid,'# Error\n');
fprintf(fid,'ERROR 0.0000001\n');
fprintf(fid,'# Collocation point depth\n');
fprintf(fid,'COLLDIST 0.0000001\n');
fprintf(fid,'# "Far field" coefficient\n');
fprintf(fid,'FARFIELD 5\n');
fprintf(fid,'# Collocation point calculation:\n#  0-approximative\n#  1-accurate\n');
fprintf(fid,'COLLCALC 0\n');
fprintf(fid,'# Interpolation method/order for velocity calculations:\n#  0-nodal\n#  1-first\n#  2-second\n');
fprintf(fid,'VELORDER 1\n\n\n');

fprintf(fid,'# RESULT REQUESTS\n#  0-no\n# 1-yes\n');
fprintf(fid,'RESULTS 1\n');
fprintf(fid,'#  1 - coefiicients\n');
fprintf(fid,'RES_COEF 1\n');
fprintf(fid,'#  2 - forces\n');
fprintf(fid,'RES_FORC 0\n');
fprintf(fid,'#  3 - geometry\n');
fprintf(fid,'RES_GEOM 0\n');
fprintf(fid,'#  4 - velocity\n');
fprintf(fid,'RES_VELO 1\n');
fprintf(fid,'#  5 - pressure\n');
fprintf(fid,'RES_PRES 1\n');
fprintf(fid,'#  6 - center points\n');
fprintf(fid,'RES_CENT 0\n');
fprintf(fid,'#  7 - dipole values\n');
fprintf(fid,'RES_DOUB 1\n');
fprintf(fid,'#  8 - source values\n');
fprintf(fid,'RES_SORC 1\n');
fprintf(fid,'#  9 - velocity components\n');
fprintf(fid,'RES_VELC 1\n');
fprintf(fid,'# 10 - mesh characteristics\n');
fprintf(fid,'RES_MESH 0\n');
fprintf(fid,'# 11 - static pressure\n');
fprintf(fid,'RES_STAT 0\n');
fprintf(fid,'# 12 - dynamic pressure\n');
fprintf(fid,'RES_DYNA 0\n');
fprintf(fid,'# 13 - manometer pressure\n');
fprintf(fid,'RES_MANO 1\n\n\n');

 
% Ler o modelo
%fid = fopen('wings_template.inp','r'); 
%string = fscanf(fid,'%c');
%fclose(fid);

%fid = fopen('wings_rewrite.inp','w');
%fprintf(fid,string);
fprintf(fid,'# GEOMETRY\n\n');
fprintf(fid,'NODES %d\n',size(NODE,1));
fprintf(fid,'%f %f %f\n',NODE');
fprintf(fid,'\nPANELS %d\n',size(PANEL,1));
fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',PANEL(1:end-sec_N+1,:)');
fprintf(fid,'%d %d %d %d %d %d %d\n',PANEL(end-sec_N+2:end,1:7)');
fclose(fid);


% Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if sim_op == 1

    system('apame_win64.exe < apame_input_v03.txt');clc

    if exist('wings_rewrite_v03.log','file')
    
        result_ID = fopen('wings_rewrite_v03.log','r');
        result = fscanf(result_ID,'%c');
        disp(result)
        fclose(result_ID);
        delete('fort.2');
        delete('wings_rewrite_v03.log');
        delete('wings_rewrite_v03.res');
    
    else
        disp('A simulação não foi bem-sucedida')
    
    end
end




