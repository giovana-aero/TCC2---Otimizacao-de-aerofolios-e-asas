% Objetivo: montar a malha de uma asa pro APAME
clc,clear,fclose('all');

% Características 
% - Asa retangular
% - Um perfil pra raiz e outro pra ponta

% Nota sobre as coordenadas de aerofólios: este script trabalha de modo que seja
% necessário que os aerofólios carregados tenham sempre o mesmo número de pontos.
% Isso não é um problema no contexto do algoritmo de otimização. No entanto, se o
% o intuito for gerar uma malha em outra aplicação, é necessário garantir esse
% requisito de números iguais de pontos. Para isso, pode-se interpolar as coordenadas
% em mãos, ou inserí-las no CST reverso e gerar o aerofólio com a quantidade desejada
% de pontos


% Dados da asa
b = 10; % Envergadura
%c_r = 0; % Corda da raiz
%c_t = 0; % Corda da ponta
%af_r = 'coordenadas.dat'; % Aerofólio da raiz (coordenadas em formato XFOIL)
%af_t = 'coordenadas_2.dat'; % Aerofólio da ponta (coordenadas em formato XFOIL)
%af_r = 'coordenadas_20_1.dat';
%af_t = 'coordenadas_20_2.dat';

% Aerofólios CST
v_ex_r = [0.04,.1,.3,.1,20,0]; v_in_r = [.01,.1,.1,.2,-10,0]; 
v_ex_t = [0.01,.1,.1,.1,10,0]; v_in_t = [.01,.1,.1,.2,10,0]; 
dat.chord = 1;
dat.BPn = 4;
dat.np = 50;
dat.N1 = 0.5;
dat.N2 = 1;
dat.p_op = 0;


% Nota: aerofólios *devem* ter o bordo de fuga fechado
%if size(af_r,1) ~= size(af_t,1),error('Número de pontos de ambos os aerofólios devem ser iguais'),end

% Dados da malha de painéis
far = b*2; % Comprimento dos painéis de trilha (a partir do bordo de fuga)
nb = 5;    % Número de seções intermediárias entre a raiz e a ponta (considerando apenas um lado da asa)

% Configurações da simulação
sim_op = 0; % Fazer a simulação no apame ao final?
v_ref = 100; % [m/s]
aoa = 0; % [graus]
rho = 1.225; % [kg/m^3]
p_atm = 101325; % [Pa]


% Carregar coordenadas dos aerofólios. Como o contorno é fechado, ignora-se o último par de coordenadas
%coo_r = dlmread(af_r); coo_r = coo_r(1:end-1,:);
%coo_t = dlmread(af_t); coo_t = coo_t(1:end-1,:);
coo_r = run_cst_TCC2(v_ex_r,v_in_r,dat); coo_r = coo_r(1:end-1,:);
coo_t = run_cst_TCC2(v_ex_t,v_in_t,dat); coo_t = coo_t(1:end-1,:);
sec_af_N = size(coo_r,1); % Número de nós por seção
sec_N = 3 + nb*2; % Número de seções transversais




% Obter nós ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NODE = zeros(sec_af_N*(1+nb),3);
NODE(1:sec_af_N,:) = [coo_t(:,1),repmat(-b/2,sec_af_N,1),coo_t(:,2)]; % Ponta esquerda da asa
%NODE = [coo_t(:,1),repmat(-b/2,sec_af_N,1),coo_t(:,2)]; % Ponta esquerda da asa

% Gerar seções intermediárias e adicionar ao lado esquerdo da asa
if nb > 0
    % Criar struct que guarda seções de asa intermediárias
    % (isto será um auxílio devido à natureza simétrica a asa)
    wing_sec.coo = []; wing_sec = repmat(wing_sec,nb,1); % Inicializar
    op_vec = linspace(1,0,2+nb); op_vec = op_vec(2:end-1); % Definição do formato da interpolação em função dos originais
    
    % Encontrar as coordenadas das seções (fazer interpolações)
    for i = 1:length(op_vec)
        wing_sec(i).coo = airfoil_interpolation(coo_r,coo_t,op_vec(i),-op_vec(i)*(b/2));
        % Adicionar ao lado esquerdo da asa
        NODE(sec_af_N*i+1:sec_af_N*(i+1),:) = wing_sec(i).coo;
    end
end

%    for i = 1:length(op_vec)
%%        wing_sec(i).coo = airfoil_interpolation(coo_r,coo_t,op_vec(i),-op_vec(i)*(b/2));
%        % Adicionar ao lado esquerdo da asa
%%        disp(NODE(sec_af_N*i+1:sec_af_N*(i+1),:))
%disp(NODE(sec_af_N*(i+1),:))
%    end

% Adicionar coordenadas do perfil da raiz
NODE = [NODE;coo_r(:,1),zeros(sec_af_N,1),coo_r(:,2)];

% Adicionar seções intermediárias ao lado direito da asa
if nb > 0
%    NODE = [NODE;zeros(sec_af_N*nb,3)];
    temp = zeros(sec_af_N*nb,3);
    k = 1;
    for i = length(op_vec):-1:1
        wing_sec(k).coo(:,2) = -wing_sec(k).coo(:,2); % Inverter o sinal da coordenada y das seções intermediárias
        temp(sec_af_N*i-sec_af_N+1:sec_af_N*i,:) = wing_sec(k).coo;
        k = k + 1;
    end
    NODE = [NODE;temp];
end

% Adicionar ponta direita
NODE = [NODE;coo_t(:,1),zeros(sec_af_N,1)+b/2,coo_t(:,2)];


% Teste: fazer um gráfico das seções da asa (raiz à ponta)
figure(2),clf,grid on,axis equal,hold on
for i = 1:(sec_N-1)
%    figure(2),clf,grid on,axis equal
    plot(NODE(sec_af_N*i-sec_af_N+1:sec_af_N*i,1),NODE(sec_af_N*i-sec_af_N+1:sec_af_N*i,3))
%    input('')
end






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
    elseif i == sec_af_N*(1+2*nb)+1 % Painel do extradorso da ponta direita(bordo de fuga)
        panel_surf(i,6:end) = [0,0,sec_af_N*(1+nb*2)+2,sec_af_N*(1+nb)+1];
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

fid = fopen('wings_rewrite.inp','w');
fprintf(fid,'APAME input file\nVERSION 3.1\n\n');
fprintf(fid,'AIRSPEED %f\n',v_ref);
fprintf(fid,'DENSITY %f\n',rho);
fprintf(fid,'PRESSURE %f\n',p_atm);
fprintf(fid,'MACH 0\n');
fprintf(fid,'CASE_NUM 1\n')
fprintf(fid,'%f\n0\n\n',aoa);
fprintf(fid,'WINGSPAN %f\n',b);
fprintf(fid,'MAC 1\n');
fprintf(fid,'SURFACE %f\n',b);
fprintf(fid,'ORIGIN *\n0 0 0\n\n');
fprintf(fid,'METHOD 0\n');
fprintf(fid,'ERROR 0.0000001\n');
fprintf(fid,'COLLDIST 0.0000001\n');
fprintf(fid,'FARFIELD 5\n');
fprintf(fid,'COLLCALC 0\n');
fprintf(fid,'VELORDER 1\n\n');
fprintf(fid,'RESULTS 1\n');
fprintf(fid,'RES_COEF 1\n');
fprintf(fid,'RES_FORC 0\n');
fprintf(fid,'RES_GEOM 0\n');
fprintf(fid,'RES_VELO 1\n');
fprintf(fid,'RES_PRES 1\n');
fprintf(fid,'RES_CENT 0\n');
fprintf(fid,'RES_DOUB 1\n');
fprintf(fid,'RES_SORC 1\n');
fprintf(fid,'RES_VELC 1\n');
fprintf(fid,'RES_MESH 0\n');
fprintf(fid,'RES_STAT 0\n');
fprintf(fid,'RES_DYNA 0\n');
fprintf(fid,'RES_MANO 1\n\n');

 
% Ler o modelo
%fid = fopen('wings_template.inp','r'); 
%string = fscanf(fid,'%c');
%fclose(fid);

%fid = fopen('wings_rewrite.inp','w');
%fprintf(fid,string);
fprintf(fid,'NODES %d\n',size(NODE,1));
fprintf(fid,'%f %f %f\n',NODE');
fprintf(fid,'\nPANELS %d\n',size(PANEL,1));
fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',PANEL(1:end-sec_N+1,:)');
fprintf(fid,'%d %d %d %d %d %d %d\n',PANEL(end-sec_N+2:end,1:7)');
fclose(fid);

if sim_op == 1
    % Fazer a simulação ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    system('apame_win64.exe < apame_input.txt');clc

    result_ID = fopen('wings_rewrite.log','r');
    result = fscanf(result_ID,'%c');
    disp(result)
    fclose(result_ID);
    delete('fort.2');
    delete('wings_rewrite.log');
    delete('wings_rewrite.res');

end
