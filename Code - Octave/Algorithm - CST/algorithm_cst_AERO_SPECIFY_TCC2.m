%% Algoritmo genético CST - Otimizador de perfis específicos

% Sobre o beta2: ele deve ter um valor que cumpra beta1+beta2>=L, onde L
% é a separação mínima entre ângulos -> beta2>=L-beta1


%% Correções a serem feitas

% - Aplicar a limitação dos ângulos em outros estágios além da geração inicial
% - Aplicar uma especificação de perfis simétricos 
% - Adicionar uma forma de mudar a precisão das variáveis https://www.mathworks.com/matlabcentral/answers/308585-how-to-generate-non-integer-random-number-in-between-two-numbers

% - Na checagem de erros: seria bom enfatizar os sinais das extensões de valores?
%   (primeira casa <= 0, segunda casa >= 0)


%%

clear,clc
fclose('all');
tic
%[y2,Fs2] = audioread('C:\Users\Guga Weffort\Documents\MATLAB\504 Finding a Treasure Box.wav');

%%

% Parâmetros do algoritmo
dat.N = 150;                              % Número de indivíduos na população (deve ser par por hora)
dat.mu = 0.05;                           % Probabilidade de mutação
dat.iter = 5;                         % Número de iterações
dat.elite = 1;                         % Aplicar elitismo?
dat.subs = 1;                          % Substituir aerofólios sem resultados? 
%dat.aero_M = zeros(dat.iter,4);  % [Tirar isto depois] (e apagar a função make_vector também)

% Dados do aerofólio a ser otimizado 
%op = 1;
dat.or_v_ex = [0.010000, 0.180472, 0.236104, 0.279989, 0.304271, 0.321324, 20.000057, 0.000000];
dat.or_v_in = [0.019999, 0.083345, 0.046642, 0.040038, 0.033295, 0.016689, -0.000479, -0.000000];

dat.symm_override = 0; % Forçar a identificação de perfil simétrico caso o algoritmo não o reconheça como tal (devido a erros de precisão)
                       % Isto efetivamente forma um perfil simétrico fazendo o espelhamento do extradorso

% Parâmetros da geometria CST
dat.BPn = length(dat.or_v_ex)-2;                  % Grau do polinômio de Bernstein (número de variáveis de design = BPn+1, desconsiderando o delta_z)
dat.np = 80;                          % Número de pontos a serem usados na geração de ordenadas
dat.p_op = 0;                % Opção da geração de pontos (1 pra cosspace, outro valor pra cosspace_half)
dat.N1 = 0.5;
dat.N2 = 1;
dat.le_R_ext1 = [-0.01,0.01,0.005]; % Limite inferior (<=0), limite superior(>=0), valor mínimo possível (pois Rle > 0)
dat.le_R_ext2 = [-0.01,0.01,0.005]; % Limite inferior (<=0), limite superior(>=0), valor mínimo possível (pois Rle > 0)
dat.A_ext1 = [-0.05,0.05]; % Limite inferior (<=0) e superior (>=0)
dat.A_ext2 = [-0.05,0.05]; % Limite inferior (<=0) e superior (>=0)
dat.B_ext1 = [-5,5]; % Limite inferior (<=0) e superior (>=0)
dat.B_ext2 = [dat.or_v_ex(end-1) + dat.or_v_in(end-1),5]; % O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
dat.chord = 1;

% Parâmetros das simulações
dat.cases = 1;                          % Número de condições de voo a serem analisadas
dat.reynolds = [1e6,1e6,1e6];           % Valores dos números de Reynolds para as simulações
dat.aoa = [0,4,8];                     % Ângulos de ataque
dat.iter_sim = [10,10,10];             % Números de iterações no XFOIL
dat.numNodes = 0;                      % Número de painéis
dat.coeff_op = ['o','!','!','o';       % Uma linha para cada condição de voo
                '!','!','!','!';
                '!','!','!','!'];
dat.coeff_val = [0.5,0.01,60,0;
                 0.5,0,0,-0.08;
                 0,0,0,-0.08];
dat.coeff_F = [1,1,1,1;
               1,1,1,1;
               1,1,1,1];
% [CL CD L/D CM] Definição da matriz dat.coeff_op
% '!' -> não usar como função objetivo
% '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD e CM)
% 'c'  -> buscar valor constante de coeficiente de momento (arbitrário)
% 'k' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val(1,4))
% 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
% em sua respectiva casa no vetor dat.coeff_val
% O vetor dat.coeff_F dá os pesos de cada função objetivo

% Checagem de erros
dat = error_check_cst_TCC2(dat);

% Template dos structs                                        
empty.v_ex = zeros(1,dat.BPn+2);
empty.v_in = zeros(1,dat.BPn+2);
empty.symm = [];
empty.aero = [];
empty.score = 0;

% Vetor que define os perfis:
% v = [ RLe A1 A2 A3 ... A(N) beta Dz]
pop = repmat(empty,dat.N,1);
chi = pop;

% Gerar população inicial
disp('<< Geração da população inicial >>')
for i = 1:dat.N
    disp(['Indivíduo ' num2str(i)])
    
    check = 0;
    while check == 0
    
        if isequal(dat.or_v_ex,dat.or_v_in) || dat.symm_override == 1 % Perfil simétrico
            % Vetor com informações do extradorso
            pop(i).v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn
                pop(i).v_ex(a) = dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex(dat.BPn+1) = dat.or_v_ex(dat.BPn+1) + rand*dat.B_ext1(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex(dat.BPn+2) = dat.or_v_ex(dat.BPn+2); %randi(dat.delta_range)*0.1*rand;   % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in = pop(i).v_ex;
        
            % Checar o raio do bordo de ataque
            if pop(i).v_ex(1) < dat.le_R_ext1(3)
                pop(i).v_ex(1) = dat.le_R_ext1(3);
                pop(i).v_in(1) = pop(i).v_ex(1);
            end
            
            % Checar separação do bordo de fuga
            if pop(i).v_ex(dat.BPn+1) + pop(i).v_in(dat.BPn+1) < dat.B_ext2(1)
                pop(i).v_ex(dat.BPn+1) = dat.B_ext2(1)/2;
                pop(i).v_in(dat.BPn+1) = dat.B_ext2(1)/2;
            end
            
            pop(i).symm = 1;
            
        else % perfil assimétrico
            % Vetor com informações do extradorso
            pop(i).v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn
                pop(i).v_ex(a) = dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex(dat.BPn+1) = dat.or_v_ex(dat.BPn+1) + rand*dat.B_ext1(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex(dat.BPn+2) = dat.or_v_ex(dat.BPn+2);  % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in(1) = dat.or_v_in(1) + rand*dat.le_R_ext2(randi(2));                     % Raio do bordo de ataque
            for a = 2:dat.BPn
                pop(i).v_in(a) = dat.or_v_in(a) + rand*dat.A_ext2(randi(2));        % Pesos intermediários
            end
            pop(i).v_in(dat.BPn+1) = (dat.B_ext2(1) - pop(i).v_ex(dat.BPn+1)) + rand*dat.B_ext2(2);     % Ângulo do bordo de fuga
            pop(i).v_in(dat.BPn+2) = dat.or_v_in(dat.BPn+2); 
            
            % Checar o raio do bordo de ataque
            if pop(i).v_ex(1) < dat.le_R_ext1(3)
                pop(i).v_ex(1) = dat.le_R_ext1(3);
            end
            if pop(i).v_in(1) < dat.le_R_ext2(3)
                pop(i).v_in(1) = dat.le_R_ext2(3);
            end
            
            % Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop(i).v_in(dat.BPn+1) < dat.B_ext2(1) - pop(i).v_ex(dat.BPn+1)
                pop(i).v_in(dat.BPn+1) = dat.B_ext2(1) - pop(i).v_ex(dat.BPn+1);
            end
            
            % Checar os pesos (soma de pesos do intradorso deve ser menor ou
            % igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop(i).v_ex(2:dat.BPn));
            sum2 = sum(pop(i).v_in(2:dat.BPn));
            if sum2 > sum1,continue,end
            
            pop(i).symm = 0;
        end
        
        % Checagem de qualidade
        check = quality(run_cst_TCC2(pop(i).v_ex,pop(i).v_in,dat),dat);
        
    end
end

if isequal(dat.or_v_ex,dat.or_v_in) || dat.symm_override == 1
    dat.symm_op = 1;
else
    dat.symm_op = 0;
end

% Gerar struct que guarda o melhor perfil de cada geração
archive = repmat(empty,dat.iter,1);

% Montar arquivo de input pro XFOIL
run_xfoil_cst_TCC2(dat,1);

%% Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for loop = 1:dat.iter
    
    % Simular os perfis e obter dados
    select = ones(1,dat.N);
    disp('<< Simulação dos aerofólios >>')
    for i = 1:dat.N
        disp(['Indivíduo ' num2str(i)])
        
        run_cst_TCC2(pop(i).v_ex,pop(i).v_in,dat,1); 
        pop(i).aero = run_xfoil_cst_TCC2(dat,2);clc % Simulação
        
        % Marcar indivíduos que não tenham convergido na simulação
        if pop(i).aero == 'n'
            select(i) = 0;
        end 
    end
        
    % Encontrar indivíduos problemáticos para já ignorá-los durante a atribuição
    % de pontuação
    % select contém os indivíduos com aero = 'n'; select2 contém os outros, que
    % convergiram nas simulações
    select = find(select == 0); select2 = 1:dat.N;
    select2(select) = [];
    
    if length(select) == dat.N
        error('Nenhum aerofólio convergiu nas simulações')
    end
    
    % Atribuir pontuações (fitnesses)
    pop = fitness_cst(pop,dat,select2);
    
    % Isto serve pra pôr todas as pontuações em um vetor
    weights = [pop.score];

    % Mostrar o melhor perfil e comparar com o original
    [~,pos] = max(weights);
    figure(1),clf
    if dat.cases == 1 % Mostrar no título o número da iteração e os dados aerodinâmicos
        plot_airfoil_cst_TCC2(2,run_cst_TCC2(pop(pos).v_ex,pop(pos).v_in,dat),loop,pop(pos)),hold on
    else % Mostrar no título apenas o número da iteração
        plot_airfoil_cst_TCC2(1,run_cst_TCC2(pop(pos).v_ex,pop(pos).v_in,dat),loop),hold on
    end
    plot_airfoil_cst_TCC2(3,run_cst_TCC2(dat.or_v_ex,dat.or_v_in,dat))
    legend('Otimizado','Original'),axis equal,grid on
    
    % Guardar o melhor perfil de cada iteração
    archive(loop) = pop(pos);
%    for i = 1:4 % [apagar isto depois]
%        % As médias contabilizam apenas os indivíduos que convergiram nas simulações
%        temp = make_vector(pop,i,select2);
%        dat.aero_M(loop,i) = mean(temp(find(temp~=0)));
%    end

    % Parar o código aqui na última iteração, já que nesse cenário o resto
    % do código é inútil
    if loop == dat.iter
        break
    end
    
    % Substituir indivíduos com pontuação nula ou negativa por aqueles com as
    % pontuações mais altas
    if dat.subs == 1
        % Agora o vetor select também inclui indivíduos que convergiram na
        % simulação, mas que são extremamente inaptos (pontuação negativa)
        select = find(weights <= 0);
        if length(select) <= length(select2) % Se o número de aerofólios com aero = 'n' for menor ou igual que o número de aerofólios com pontuação
            % Preencher o vetor ind com os índices dos indivíduos de maior pontuação
            ind = zeros(1,length(select));
            temp = weights;
            for i = 1:length(ind)
                [~,ind(i)] = max(temp);
                temp(ind(i)) = -Inf;
            end
            % Substituir os indivíduos
            k = 1;
            for i = [select]
                pop(i) = pop(ind(k));
                k = k + 1;
            end
            
        else 
            % Preencher o vetor ind com os índices dos indivíduos de maior pontuação (que acabam sendo todos aqueles apontados por select2)
            ind = select2;
            % Substituir os indivíduos
            k = 1;
            for i = [select]
                pop(i) = pop(ind(k));
                k = k + 1;
                if k > length(ind),k = 1;end
            end
        end
    end
    
    % Escolher membros da população pra reprodução
    c = 1;
    disp('<< Reprodução >>')
    for f = 1:dat.N/2
        fprintf('%.2f%% completo\n',f/(dat.N/2)*100)
        
        % Crossover, gerando dois filhos de cada dois pais
        % Aqui ignora-se os últimos dois parâmetros (f e delta_z) porque eles
        % sempre serão iguais
        check = 0;
        while check == 0
            
            % Isto seleciona dois pais por meio de uma seleção via roleta
            % (indivíduos com pesos maiores têm  mais chance de serem selecionados)
            par = [0,0]; % Vetor que indica a numeração dos pais escolhidos
            par(1) = selection_crossover(weights);
            par(2) = selection_crossover(weights);
            
            % v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
            if pop(par(1)).symm == 0 && pop(par(2)).symm == 0 % Se ambos forem assimétricos
                n = randi([1,4]);
                switch n
                    case 1 % Trocar os extradorsos e intradorsos inteiros
                        chi(c).v_ex = pop(par(1)).v_ex;
                        chi(c).v_in = pop(par(2)).v_in;
                        chi(c+1).v_ex = pop(par(2)).v_ex;
                        chi(c+1).v_in = pop(par(1)).v_in;
                        
                    case 2 % Trocar o raio do bordo de ataque
                        chi(c).v_ex = [pop(par(1)).v_ex(1),pop(par(2)).v_ex(2:end)];
                        chi(c).v_in = [pop(par(1)).v_in(1),pop(par(2)).v_in(2:end)];
                        chi(c+1).v_ex = [pop(par(2)).v_ex(1),pop(par(1)).v_ex(2:end)];
                        chi(c+1).v_in = [pop(par(2)).v_in(1),pop(par(1)).v_in(2:end)];
                        
                    case 3 % Trocar os pesos intermediários
                        if dat.BPn == 2
                            op = 1
                        else
                            op = randi([1 2]);
                        end
                        if op == 1 % Trocar tudo    
                            chi(c).v_ex = [pop(par(2)).v_ex(1),pop(par(1)).v_ex(2:dat.BPn),pop(par(2)).v_ex(dat.BPn+1:end)];
                            chi(c).v_in = [pop(par(2)).v_in(1),pop(par(1)).v_in(2:dat.BPn),pop(par(2)).v_in(dat.BPn+1:end)];
                            chi(c+1).v_ex = [pop(par(1)).v_ex(1),pop(par(2)).v_ex(2:dat.BPn),pop(par(1)).v_ex(dat.BPn+1:end)];
                            chi(c+1).v_in = [pop(par(1)).v_in(1),pop(par(2)).v_in(2:dat.BPn),pop(par(1)).v_in(dat.BPn+1:end)];
                        else % Trocar cortes
                            num1 = randi(2:dat.BPn);
                            num2 = randi(2:dat.BPn);
                            temp1_1 = [pop(par(1)).v_ex(2:num1),pop(par(2)).v_ex(num1+1:dat.BPn)];
                            temp1_2 = [pop(par(1)).v_in(2:num2),pop(par(2)).v_in(num2+1:dat.BPn)];
                            temp2_1 = [pop(par(2)).v_ex(2:num1),pop(par(1)).v_ex(num1+1:dat.BPn)];
                            temp2_2 = [pop(par(2)).v_in(2:num2),pop(par(1)).v_in(num2+1:dat.BPn)];
                            chi(c).v_ex = [pop(par(1)).v_ex(1),temp2_1,pop(par(1)).v_ex(dat.BPn+1:end)];
                            chi(c).v_in = [pop(par(1)).v_in(1),temp2_2,pop(par(1)).v_in(dat.BPn+1:end)];
                            chi(c+1).v_ex = [pop(par(2)).v_ex(1),temp1_1,pop(par(2)).v_ex(dat.BPn+1:end)];
                            chi(c+1).v_in = [pop(par(2)).v_in(1),temp1_2,pop(par(2)).v_in(dat.BPn+1:end)]; 
                        end
                        
                    case 4 % Trocar os ângulos do bordo de fuga
                        chi(c).v_ex = [pop(par(2)).v_ex(1:dat.BPn),pop(par(1)).v_ex(dat.BPn+1),pop(par(2)).v_ex(dat.BPn+2)];
                        chi(c).v_in = [pop(par(2)).v_in(1:dat.BPn),pop(par(1)).v_in(dat.BPn+1),pop(par(2)).v_in(dat.BPn+2)];
                        chi(c+1).v_ex = [pop(par(1)).v_ex(1:dat.BPn),pop(par(2)).v_ex(dat.BPn+1),pop(par(1)).v_ex(dat.BPn+2)];
                        chi(c+1).v_in = [pop(par(1)).v_in(1:dat.BPn),pop(par(2)).v_in(dat.BPn+1),pop(par(1)).v_in(dat.BPn+2)];
                end
                
                if n == 1
                    % Checar a separação dos bordos de fuga. Se não cumprirem o 
                    % requisito de separação, alterar o ângulo do bordo de fuga
                    % do intradorso
                    if chi(c).v_in(dat.BPn+1) < (dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1))
                        chi(c).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1);
                    end
                    if chi(c+1).v_in(dat.BPn+1) < (dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn+1))
                        chi(c+1).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn+1);
                    end
                    % Decisão de alterar o intradorso em base de uma nota na página
                    % 57(87) do Raymer (2018)
                end
                
                % Checar os pesos
                sum1 = sum(chi(c).v_ex(2:dat.BPn));
                sum2 = sum(chi(c).v_in(2:dat.BPn));
                if sum2 > sum1,continue,end
                sum1 = sum(chi(c+1).v_ex(2:dat.BPn));
                sum2 = sum(chi(c+1).v_in(2:dat.BPn));
                if sum2 > sum1,continue,end
                
                % Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                % for usado como variável)
                %chi(c).v_in(dat.BPn) = -chi(c).v_ex(dat.BPn)*chi(c).v_ex(dat.BPn-1)/chi(c).v_in(dat.BPn-1);
                %chi(c+1).v_in(dat.BPn) = -chi(c+1).v_ex(dat.BPn)*chi(c+1).v_ex(dat.BPn-1)/chi(c+1).v_in(dat.BPn-1);
            
                chi(c).symm = 0;
                chi(c+1).symm = 0;
        
            elseif pop(par(1)).symm == 1 && pop(par(2)).symm == 1 % Se ambos forem simétricos
				n = randi([1,3]);
                switch n                        
                    case 1 % Trocar o raio do bordo de ataque
                        chi(c).v_ex = [pop(par(1)).v_ex(1),pop(par(2)).v_ex(2:end)];
                        chi(c).v_in = chi(c).v_ex;
                        chi(c+1).v_ex = [pop(par(2)).v_ex(1),pop(par(1)).v_ex(2:end)];
                        chi(c+1).v_in = chi(c+1).v_ex;
                        
                    case 2 % Trocar os pesos intermediários
                        if dat.BPn == 2
                            op = 1
                        else
                            op = randi([1 2]);
                        end
                        if op == 1 % Trocar tudo    
                            chi(c).v_ex = [pop(par(2)).v_ex(1),pop(par(1)).v_ex(2:dat.BPn),pop(par(2)).v_ex(dat.BPn+1:end)];
                            chi(c).v_in = chi(c).v_ex;
                            chi(c+1).v_ex = [pop(par(1)).v_ex(1),pop(par(2)).v_ex(2:dat.BPn),pop(par(1)).v_ex(dat.BPn+1:end)];
                            chi(c+1).v_in = chi(c+1).v_ex;
                        else % Trocar cortes
                            num1 = randi(2:dat.BPn);
                            temp1_1 = [pop(par(1)).v_ex(2:num1),pop(par(2)).v_ex(num1+1:dat.BPn)];
                            temp2_1 = [pop(par(2)).v_ex(2:num1),pop(par(1)).v_ex(num1+1:dat.BPn)];
                            chi(c).v_ex = [pop(par(1)).v_ex(1),temp2_1,pop(par(1)).v_ex(dat.BPn+1:end)];
                            chi(c).v_in = chi(c).v_ex;
                            chi(c+1).v_ex = [pop(par(2)).v_ex(1),temp1_1,pop(par(2)).v_ex(dat.BPn+1:end)];
                            chi(c+1).v_in = chi(c+1).v_ex;
                        end
                        
                    case 3 % Trocar os ângulos do bordo de fuga
                        chi(c).v_ex = [pop(par(2)).v_ex(1:dat.BPn),pop(par(1)).v_ex(dat.BPn+1),pop(par(2)).v_ex(dat.BPn+2)];
                        chi(c).v_in = [pop(par(2)).v_in(1:dat.BPn),pop(par(1)).v_in(dat.BPn+1),pop(par(2)).v_in(dat.BPn+2)];
                        chi(c+1).v_ex = [pop(par(1)).v_ex(1:dat.BPn),pop(par(2)).v_ex(dat.BPn+1),pop(par(1)).v_ex(dat.BPn+2)];
                        chi(c+1).v_in = [pop(par(1)).v_in(1:dat.BPn),pop(par(2)).v_in(dat.BPn+1),pop(par(1)).v_in(dat.BPn+2)];
                end
                
                chi(c).symm = 1;
                chi(c+1).symm = 1;
            
            else % Se um for simétrico e o outro for assimétrico 
                % Transformar o simétrico em um assimétrico e vice-versa
                % chi(c) é simétrico e chi(c+1) é assimétrico
                if pop(par(1)).symm == 0 
                    chi(c).v_ex = pop(par(1)).v_ex;
                    chi(c).v_in = pop(par(1)).v_ex;
                    chi(c+1).v_ex = pop(par(2)).v_ex;
                    chi(c+1).v_in = pop(par(1)).v_in;
                else
                    chi(c).v_ex = pop(par(2)).v_ex;
                    chi(c).v_in = pop(par(2)).v_ex;
                    chi(c+1).v_ex = pop(par(1)).v_ex;
                    chi(c+1).v_in = pop(par(2)).v_in;
                end
                         
                
                % Checar a separação dos bordos de fuga. Se não cumprirem o 
                % requisito de separação, alterar o ângulo do bordo de fuga
                % do intradorso
                if chi(c).v_in(dat.BPn+1) < (dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1))
                    chi(c).v_ex = dat.B_ext2(1)/2;
                    chi(c).v_in = dat.B_ext2(1)/2;
                end
                if chi(c+1).v_in(dat.BPn+1) < (dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn+1))
                    chi(c+1).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn+1);
                end
                % Decisão de alterar o intradorso em base de uma nota na página
                % 57(87) do Raymer (2018)
                
                
                % Checar os pesos
%                sum1 = sum(chi(c).v_ex(2:dat.BPn));
%                sum2 = sum(chi(c).v_in(2:dat.BPn));
%                if sum2 > sum1,continue,end
                sum1 = sum(chi(c+1).v_ex(2:dat.BPn));
                sum2 = sum(chi(c+1).v_in(2:dat.BPn));
                if sum2 > sum1,continue,end
                
                chi(c).symm = 1;
                chi(c+1).symm = 0;
                
            end
            
			% Checagem de qualidade
            check = quality(run_cst_TCC2(chi(c).v_ex,chi(c).v_in,dat),dat);
            if check == 0,continue,end
            check = quality(run_cst_TCC2(chi(c+1).v_ex,chi(c+1).v_in,dat),dat); 
			
        end
        c = c + 2;
    end
    
    % Mutação
    select1 = (rand(size(pop)) <= dat.mu);
    if sum(select1) ~= 0
		clc,disp('<< Mutação >>')
        select2 = find(select1 == 1);
        k = 1;
        for i = select2
            s = select2(k);
            disp(['Indivíduo ' num2str(k) ' de ' num2str(length(select2))])
            
            
            
%            pop(i).v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));    % Raio do bordo de ataque
%            for a = 2:dat.BPn
%                pop(i).v_ex(a) = dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));        % Pesos intermediários
%            end
%            pop(i).v_ex(dat.BPn+1) = dat.or_v_ex(dat.BPn+1) + rand*dat.B_ext1(randi(2));     % Ângulo do bordo de fuga
%            pop(i).v_ex(dat.BPn+2) = dat.or_v_ex(dat.BPn+2);  % delta_z
            
            
            
            check = 0;
            while check == 0
                
                temp = chi(s);
                if temp.symm == 1 % Caso o perfil seja simétrico
                    n = [1,4,5](randi(3));
                else % Caso o perfil seja assimétrico
                    n = randi(5);
                end
                switch n
                    case 1 % Alterar o raio do bordo de ataque
                        if temp.symm == 1
                            P = 1;
                        else
                            p = randi([1,4]);
                        end
                        if p == 1 % Mudar ambos para o mesmo valor
                            temp.v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));
                            temp.v_in(1) = temp.v_ex(1);
                        elseif p == 2 % Mudar ambos independentemente 
                            temp.v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));
                            temp.v_in(1) = dat.or_v_in(1) + rand*dat.le_R_ext2(randi(2));
                        elseif p == 3 % Mudar do extradorso
                            temp.v_ex(1) = dat.or_v_ex(1) + rand*dat.le_R_ext1(randi(2));
                        else % Mudar do intradorso
                            temp.v_in(1) = dat.or_v_in(1) + rand*dat.le_R_ext2(randi(2));
                        end
                        
                        % Checar o raio do bordo de ataque (simétricos)
                        if temp.symm == 1 && temp.v_ex(1) < dat.le_R_ext1(3)
                            temp.v_ex(1) = dat.le_R_ext1(3);
                            temp.v_in(1) = temp.v_ex(1);
                        end
                        
                        % Checar o raio do bordo de ataque (assimétricos)
                        if temp.symm == 0 && temp.v_ex(1) < dat.le_R_ext1(3)
                            temp.v_ex(1) = dat.le_R_ext1(3);
                        end
                        if temp.symm == 0 && temp.v_in(1) < dat.le_R_ext2(3)
                            temp.v_in(1) = dat.le_R_ext2(3);
                        end
                        
                    case 2 % Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
%                        num = rand(1,dat.BPn-1)*0.1;
                        for a = 2:dat.BPn
                            temp.v_ex(a) = dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));
                        end
                        
                    case 3 % Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
%                        num = rand(1,dat.BPn-1)*0.1;
                        for a = 2:dat.BPn
                            temp.v_in(a) = dat.or_v_in(a) + rand*dat.A_ext2(randi(2));
                        end
                        
                    case 4 % Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                        if temp.symm == 1
%                            num = rand(1,dat.BPn-1)*0.1;
                            for a = 2:dat.BPn
                                temp.v_ex(a) =dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));
                                temp.v_in(a) = temp.v_ex(a);
                            end
                        else
%                            num = rand(1,dat.BPn-1)*0.1;
                            for a = 2:dat.BPn
                                temp.v_ex(a) = dat.or_v_ex(a) + rand*dat.A_ext1(randi(2));
                            end
                            
%                            num = rand(1,dat.BPn-1)*0.1;
                            for a = 2:dat.BPn
                                temp.v_in(a) = dat.or_v_in(a) + rand*dat.A_ext2(randi(2));
                            end
                        end
                        
                    case 5 % Alterar o ângulo do bordo de fuga 
                        if temp.symm == 1 % Se for simétrico
                            temp.v_ex(dat.BPn+1) = dat.or_v_ex(dat.BPn+1) + rand*dat.B_ext1(randi(2));
                            temp.v_in(dat.BPn+1) = temp.v_ex(dat.BPn+1);
                            
                            % Checar separação do bordo de fuga
                            if temp.v_ex(dat.BPn+1) + temp.v_in(dat.BPn+1) < dat.B_ext2(1)
                                temp.v_ex(dat.BPn+1) = dat.B_ext2(1)/2;
                                temp.v_in(dat.BPn+1) = dat.B_ext2(1)/2;
                            end
                        else % Se for assimétrico
                            temp.v_ex(dat.BPn+1) = dat.or_v_ex(dat.BPn+1) + rand*dat.B_ext1(randi(2));
                            temp.v_in(dat.BPn+1) = (dat.B_ext2(1) - temp.v_ex(dat.BPn+1)) + rand*dat.B_ext2(2)
                            
                            % Checar separação do bordo de fuga (beta2>=L-beta1)
                            if temp.v_in(dat.BPn+1) < dat.B_ext2(1) - temp.v_ex(dat.BPn+1)
                                temp.v_in(dat.BPn+1) = dat.B_ext2(1) - temp.v_ex(dat.BPn+1);
                            end
                        end
                        
                end
                 
                % Checar os pesos (soma de pesos do intradorso deve ser menor ou
                % igual à soma de pesos do extradorso)
                sum1 = sum(temp.v_ex(2:dat.BPn));
                sum2 = sum(temp.v_in(2:dat.BPn));
                if sum2 > sum1
                    continue
                end
                
                % Checagem de qualidade
                check = quality(run_cst_TCC2(temp.v_ex,temp.v_in,dat),dat); 
                


                
                % Debugging: checar se o requisito de valor mínimo dos raios do bordo de ataque está sendo cumprido
                if temp.v_ex(1) < dat.le_R_ext1(3) || temp.v_in(1) < dat.le_R_ext2(3)
                    error('Valor mínimo de raio de bordo de ataque não cumprido')
                end
                
            end
            chi(s) = temp;
            k = k + 1;
        end
    end

    % substituir a população inicial pelos filhos
    pop = chi;
    
    % Aplicar elitismo
    if dat.elite == 1
        % Passar o melhor indivíduo pra nova população
        pop(1) = archive(loop);
        pop(1).score = 0;
    end
    
end



%% Final ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Fazer simulação do aerofólio original para fins de comparação
run_cst_TCC2(dat.or_v_ex,dat.or_v_in,dat,1);
original.aero = run_xfoil_cst_TCC2(dat,2);

clc,    %sound(y2,Fs2)
t = toc;
min = fix(t/60); s = rem(t,60);
fprintf('Tempo: %d min e %.2f s\n', min, s)

% Fazer gráficos dos coeficientes dos melhores indivíduos de cada geração
for i = 1:dat.cases
    aero_m = zeros(dat.iter,4);
    for j = 1:dat.iter
        aero_m(j,:) = archive(j).aero(i,:);
    end
    figure(i+1),clf,hold on,grid on
    h1 = plot(1:dat.iter,aero_m(:,1)','g-*');
    h2 = plot(1:dat.iter,aero_m(:,2)','r-*');
    [hax,h4,h3] = plotyy(1:dat.iter,aero_m(:,4)',1:dat.iter,aero_m(:,3)');
    set(h3,'color','k'),set(h3,'marker','*'),set(h4,'color','b'),set(h4,'marker','*')
    lines = [h1,h2,h3,h4]; legend(lines,'CL','CD','CM','L/D');
    legend([h1,h2,h3,h4],'CL','CD','L/D','CM');
    xlabel(hax(1),'Iteração')
    ylabel(hax(1),'CL, CD e CM')
    ylabel(hax(2),'L/D')
    title(['Melhores resultados - Condição de voo ' num2str(i) ':  Re ' num2str(dat.reynolds(i)) ', AoA ' num2str(dat.aoa(i)) '°'])
    
    % Trocar espaçador decimal
    xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
    new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
    set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

end

%% Pegar o struct de arquivo e imprimir todos
fprintf('Grau do polinômio: %d\n',dat.BPn)
for i = 1:length(archive)
    disp(['<< Iteração ' num2str(i) ' >>'])
    
    fprintf('v_ex = [%.4f, ', archive(i).v_ex(1))
    for j = 2:(length(pop(1).v_ex)-3)
        fprintf('%.4f, ',archive(i).v_ex(j))
    end
    fprintf('%.4f, ', archive(i).v_ex(end-2))
    fprintf('%.4f, ', archive(i).v_ex(end-1))
    fprintf('%.4f];\n', archive(i).v_ex(end))
    
    fprintf('v_in = [%.4f, ', archive(i).v_in(1))
    for j = 2:(length(pop(1).v_in)-3)
        fprintf('%.4f, ',archive(i).v_in(j))
    end
    fprintf('%.4f, ', archive(i).v_in(end-2))
    fprintf('%.4f, ', archive(i).v_in(end-1))
    fprintf('%.4f];\n', archive(i).v_in(end))
    for j = 1:dat.cases
        fprintf('Condição de voo %d (CL,CD,L/D,CM)\nOriginal: ',j),disp(original.aero(j,:))
        fprintf('Otimizado: '),disp(archive(i).aero(j,:))
    end
    fprintf('Pontuação: '),disp(archive(i).score)
    disp('')
end




% Fazer gráfico das médias dos coeficientes [apagar isto depois]
%figure(10),clf
%h1 = plot(1:dat.iter,dat.aero_M(:,1)','g-*');hold on,grid on
%h2 = plot(1:dat.iter,dat.aero_M(:,2)','r-*');
%[hax,h4,h3] = plotyy(1:dat.iter,dat.aero_M(:,4)',1:dat.iter,dat.aero_M(:,3)');
%set(h3,'color','k'),set(h3,'marker','*'),set(h4,'color','b'),set(h4,'marker','*')
%legend([h1,h2,h3,h4],'CL','CD','L/D','CM');
%xlabel(hax(1),'Iteração')
%ylabel(hax(1),'CL, CD e CM')
%ylabel(hax(2),'L/D')
%title('Médias dos coeficientes em cada iteração')

% Template de comando pra pegar coordenadas e um gráfico de determinado perfil
%clf,plot_airfoil_cst_TCC2(4,run_cst_TCC2(v_ex,v_in,dat,0)),grid on,axis equal

%figure(2),saveas(gcf,'fig1.png')
%figure(3),saveas(gcf,'fig2.png')
%dat.or_v_ex = [0.019257, 0.186300, 0.249870, 0.157309, 0.235440, 0.188924, 12.410294, -0.000000]; % Este é um perfil NACA 4 dígitos
%dat.or_v_in = [0.011473, 0.115911, 0.093905, 0.101085, 0.080645, 0.078528, 4.450847, -0.000000];
%figure(10),saveas(gcf,'aoba.png')