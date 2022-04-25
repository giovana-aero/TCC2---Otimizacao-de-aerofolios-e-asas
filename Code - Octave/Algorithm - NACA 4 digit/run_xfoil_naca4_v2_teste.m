% Modificação (TCC2): capacidade para simular múltiplas condições de voo
clc,clear

% Parâmetros do algoritmo
dat.N = 500;                              % Número de indivíduos na população (deve ser par por hora)
dat.mu = 0.05;                           % Probabilidade de mutação
dat.iter = 5;                           % Número de iterações
dat.aero = zeros(dat.iter,4);          % Matriz usada no gráfico final
dat.aero_m = dat.aero;
dat.elite = 1;                         % Aplicar elitismo?
dat.subs = 1;                          % Substituir aerofólios sem resultados? 

% Parâmetros da geometria
dat.m_ext = [0 9];
dat.p_ext = [1 9];
dat.t_ext = [10 30];

% Parâmetros das simulações
dat.cases = 2;                          % Número de condições de voo a serem analisadas
dat.reynolds = [1e6,1e6,1e7];           % Valores dos números de Reynolds para as simulações
dat.aoa = [0,2,0];                     % Ângulos de ataque
dat.iter_sim = [50,50,50];             % Números de iterações no XFOIL
%dat.coeff_op = ['!','!','!','^';       % Uma linha para cada condição de voo
%                '!','!','!','^';
%                '!','!','!','^'];
%dat.coeff_val = [0.6,0,120,0;
%                 0,0,0,0;
%                 0,0,0,0];
%dat.coeff_F = [1,1,1,1;
%               1,1,1,1;
%               1,1,1,1];

naca_num = '0013';



reynolds = dat.reynolds;
aoa = dat.aoa;
iter_sim = dat.iter_sim;
%coeff_op = dat.coeff_op;
cases = dat.cases;


%fclose('all');





% Apagar arquivos caso existam
for i = 1:cases
    a_polar = ['polar' num2str(i) '.txt'];
    if (exist(a_polar,'file'))
        delete(a_polar);
    end
end

%if (exist('xfoil_input.txt','file'))
%    delete('xfoil_input.txt');
%end


%% Criação do arquivo de input do Xfoil
fid = fopen('xfoil_input.txt','w');


% Mudar uma opção dos gráficos do XFOIL (desativar a aparição 
% da janela com o desenho da simulação)
fprintf(fid,'PLOP\n');
fprintf(fid,'G\n\n');

% Gerar o aerofólio
fprintf(fid, 'NACA %s\n', naca_num);

% Simulações 
fprintf(fid,'OPER\n');
if cases == 1 % Apenas uma condição de voo
    
    fprintf(fid,'VISC %d\n', reynolds(1)); % Aplicar o modo viscoso e determinar número de Reynolds 
    fprintf(fid,'PACC\n');              % Estabelecer arquivo de output
    fprintf(fid,['polar' num2str(1) '.txt\n\n']);
    fprintf(fid,'ITER %d\n', iter_sim(1)); % Mudar o número de iterações      
    fprintf(fid,'ALFA %f\n', aoa(1));        % Estabelecer ângulo de ataque

else % Duas ou mais condições de voo

    for i = 1:cases
        
        if i == 1
            fprintf(fid,'VISC %d\n', reynolds(i)); % Aplicar o modo viscoso e determinar número de Reynolds
        elseif i ~= 1 && reynolds(i) ~= reynolds(i-1)
            fprintf(fid,'RE %d\n',reynolds(i)); % Mudar número de Reynolds (após primeira simulação)
%        else
%            fprintf(fid,'RE %d\n',reynolds(i)); % Mudar número de Reynolds (após primeira simulação)
        end
        fprintf(fid,'PACC\n');              % Estabelecer arquivo de output
        fprintf(fid,['polar' num2str(i) '.txt\n\n']);
        if i == 1
            fprintf(fid,'ITER %d\n', iter_sim(i)); % Mudar o número de iterações    
        elseif i~= 1 && iter_sim(i) ~= iter_sim(i-1)
            fprintf(fid,'ITER %d\n', iter_sim(i));
        end
        fprintf(fid,'ALFA %f\n', aoa(i));        % Estabelecer ângulo de ataque
        fprintf(fid,'PACC\n'); % Fechar a polar para começar a próxima simulação
        
    end
end

% Fechar arquivo
fprintf(fid,'\nQUIT\n');
fclose(fid);

% Executar XFOIL com o arquivo de entrada
cmd = 'xfoil.exe < xfoil_input.txt' ;
system(cmd);


%% Ler o arquivo de output
aero = zeros(cases,4);
for i = 1:cases

    % Contar o número de linhas
    fidpolar = fopen(['polar' num2str(i) '.txt']); 
    tline = fgetl(fidpolar);
    nl = 0;
    while ischar(tline)	
        nl = nl+1; % Número de linhas
        tline = fgetl(fidpolar);
    end
    fclose(fidpolar);
    
    if nl == 13
        fidpolar = fopen(['polar' num2str(i) '.txt']);
        dataBuffer = textscan(fidpolar,'%f %f %f %f %f %f %f','HeaderLines',12,...
            'CollectOutput',1,...
            'Delimiter','');
        fclose(fidpolar);
        
        % Valores dos coeficientes
        CL = dataBuffer{1,1}(1,2);
        CD = dataBuffer{1,1}(1,3);
        CM = dataBuffer{1,1}(1,5);
%        delete(a_polar);
        
        aero(i,:) = [CL CD CL/CD CM];
        
    else
        aero = 'n'; break
    end

    
end

disp(''),disp(aero)

%fclose(fidpolar);

