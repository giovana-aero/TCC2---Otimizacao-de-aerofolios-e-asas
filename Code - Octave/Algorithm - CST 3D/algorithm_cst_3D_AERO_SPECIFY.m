%% Algoritmo genético
% Otimização de asas trapezoidais simples ou duplas, com perfis CST, feita
% em termos de uma geometria inicial específica


% Nota pra posterioridade: grandes quantidades de pontos nos aerofólios (np)
% tendem a fazer o apame travar de vez em quando



% A resolver:
% - na extensão de b1, inserir mais uma casa pro valor mínimo possível
% - Checar se a opção 'L' da corda do meio está sendo satisfeita nos crossovers
%   e nas mutações

% - Mutação de aerofólios: mutação dos ângulos de bordo de fuga estão erradas (case 5)



clear,clc,close('all');
fclose('all');
tic


% Parâmetros do algoritmo
dat.N = 450;                              % Número de indivíduos na população
dat.mu = 0.05;                           % Probabilidade de mutação (definida entre zero e um)
dat.iter = 5;                           % Número de iterações
dat.elite = 1;                         % Aplicar elitismo?
dat.subs = 1;                          % Substituir asas sem resultados? (ver ainda se isto será necessário)
dat.aero_M = zeros(dat.iter,4);  % [Tirar isto depois] (e apagar a função make_vector também)

% Dados da asa a ser otimizada
% Planta
dat.half = 0;
dat.type = 1; % (0 -> trapezoidal simples, 1 -> trapezoidal dupla) 
dat.or_b = 14;
dat.or_b1 = 12;
dat.or_c_r = 1.5;
dat.or_c_m = 1;
dat.or_c_t = 0.5;
dat.or_sweep = 0;
dat.or_sweep1 = 2;
dat.or_sweep2 = 20;
% Aerofólio da raiz
dat.or_v_ex_r = [0.0150, 0.2500, 0.2750, 0.1750, 8.0000, 0.0000];
dat.or_v_in_r = [0.0250, 0.1000, 0.1000, 0.1000, 2.0000, 0.0000];
dat.symm_override_r = 0;
% Aerofólio do meio
dat.or_v_ex_m = [0.0200, 0.1500, 0.1250, 0.1500, 8.0000, 0.0000];
dat.or_v_in_m = [0.0050, 0.0000, 0.0500, 0.0000, 2.0000, 0.0000];
dat.symm_override_m = 0;
% Aerofólio da ponta
dat.or_v_ex_t = [0.0100, 0.2000, 0.2000, 0.0000, 10.0000, 0.0000];
dat.or_v_in_t = [0.0050, 0.0000, 0.0750, 0.0000, 5.0000, 0.0000];
dat.symm_override_t = 0;
dat.or_tw_t = 0;

% Parâmetros da geometria: planta da asa
%dat.planf_op = 0.5; % Proporção de asas trapezoidais simples e bitrapezoidais (0->todas trapezoidais simples, 1->todas bitrapezoidais)
dat.b_ext = [-0,0,1]; % Envergadura completa [m] (limite inferior, limite superior, valor mínimo possível)
dat.b1_ext = [-0.5,0.5,1,0.5]; % Envergadura da raiz ao meio [m] (asas bitrapezoidais apenas) (valor mínimo,valor máximo,valor mínimo possível,separação mínima da ponta da asa (considerando apenas uma metade)) 
dat.c_r_ext = [-0.1,0.1,0.5]; % Corda da raiz [m] (limite inferior, limite superior, valor mínimo possível)
dat.c_m_ext = [-0.1,0.1,0.5]; % Corda do meio [m] (asas bitrapezoidais apenas) (limite inferior, limite superior, valor mínimo possível) (a opção 'L' força o formato trapezoidal simples)
dat.c_t_ext = [-0.1,0.1,0.1]; % Corda da ponta [m] (limite inferior, limite superior, valor mínimo possível)
dat.sweep_ext = [-0,0]; % Enflechamento de asas trapezoidais simples (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep1_ext = [-5,5]; % Enflechamento da primeira seção de asas trapezoidais duplas (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
dat.sweep2_ext = [-5,5]; % Enflechamento da segunda seção de asas trapezoidais duplas (opção 'Z' faz com que a linha c/2 tenha enflechamento zero)
% Parâmetros da geometria: aerofólio da raiz
dat.BPn_r = length(dat.or_v_ex_r)-2;                  % Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_r = 0.5;
dat.N2_r = 1;
dat.le_R_ext1_r = [-0.0,0.0,0.005]; % Limite inferior, limite superior, valor mínimo possível (pois Rle > 0)
dat.le_R_ext2_r = [-0.0,0.0,0.005];
dat.A_ext1_r = [-0.,0.];
dat.A_ext2_r = [-0.,0.];
dat.B_ext1_r = [0,0]; % Limites inferior e superior
dat.B_ext2_r = [dat.or_v_ex_r(end-1) + dat.or_v_in_r(end-1),0]; % O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
% Parâmetros da geometria: aerofólio do meio (asas bitrapezoidais apenas) 
dat.BPn_m = length(dat.or_v_ex_m)-2;                  % Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_m = 0.5;
dat.N2_m = 1;
dat.le_R_ext1_m = [-0.0,0.0,0.005]; % (admite opção 'L')
dat.le_R_ext2_m = [-0.0,0.0,0.005];
dat.A_ext1_m = [-0.,0.];
dat.A_ext2_m = [-0.,0.];
dat.B_ext1_m = [-0,0]; % Limites inferior e superior
dat.B_ext2_m = [dat.or_v_ex_m(end-1) + dat.or_v_in_m(end-1),0]; % O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
% Parâmetros da geometria: aerofólio da ponta
dat.BPn_t = length(dat.or_v_ex_t)-2;                  % Grau do polinômio de Bernstein (número de variáveis de design = BPn+1 mais o delta_z)
dat.N1_t = 0.5;
dat.N2_t = 1;
dat.le_R_ext1_t = [-0.0,0.0,0.005]; 
dat.le_R_ext2_t = [-0.0,0.0,0.005];
dat.A_ext1_t = [-0.,0.];
dat.A_ext2_t = [-0.,0.];
dat.B_ext1_t = [-0,0]; % Limites inferior e superior
dat.B_ext2_t = [dat.or_v_ex_t(end-1) + dat.or_v_in_t(end-1),0]; % O primeiro número é a separação mínima do extradorso, o segundo é o limite superior
dat.tw_t_ext = [0,0];

% Parâmetros da malha
dat.np = 30; % Número de pontos na geração de ordenadas nos aerofólios
dat.np_op = 1; % 1 -> cosspace, 0 -> cosspace_half
dat.nb = [2,1]; % Número de seções intermediárias (raiz/ponta) [número de seções,0] ou [concentração por metro,1]
dat.nb1 = 'L'; % Número de seções intermediárias (raiz/meio) (asas bitrapezoidais apenas) (opção 'L' faz com que nb1 e nb2 sejam uniformemente determinados ao longo da envergadura)
dat.nb2 = []; % Número de seções intermediárias (meio/ponta) (asas bitrapezoidais apenas)

% Parâmetros das simulações
dat.cases = 1;                          % Número de condições de voo a serem analisadas
dat.v_ref = [100,100,100]; % Velocidades de referência [m/s] 
dat.rho = [1.225,1.225,1.225]; % Densidades do ar [kg/m^3] 
dat.p_atm = [101325,101325,101325]; % Pressões do ar [Pa] (irrelevante neste algoritmo)
dat.mach = [0,0,0.]; % Números de Mach
dat.reynolds = [1e6,1e6,1e6];           % Valores dos números de Reynolds para as simulações (irrelevante neste algoritmo))
dat.aoa = [5,0,4];                     % Ângulos de ataque
dat.coeff_op = ['!','^','!','!';       % Uma linha para cada condição de voo
                '!','!','!','!';
                '!','!','!','!'];
dat.coeff_val = [0.3,7e-3,90,-1e-1;
                 0.5,0,0,-0.08;
                 0,0,0,-0.08];
dat.coeff_F = [1,1,1,1;
               1,1,1,1;
               1,1,1,1];
% [CL CD L/D CM] Definição de cada linha da matriz dat.coeff_op
% '!' -> não usar como função objetiva
% '^' -> procurar por um valor máximo (CL e L/D) ou valor mínimo (CD)
% 'c' -> buscar valor constante de coeficiente de momento (arbitrário)
% 'k' -> buscar valor constante de coeficiente de momento (específico, de dat.coeff_val(1,4))
% 'o' -> procurar por um valor específico (qualquer um dos parâmetros). Nesse caso, definir o valor
% em sua respectiva casa na matriz dat.coeff_val
% 'q' -> procurar por um valor específico de força de sustentação, força de arrasto
% ou momento de arfagem (CL, CD e CM). Nesse caso, definir o valor em sua casa na matriz dat.coeff_val
% '#' -> procurar por um valor máximo (L) ou mínimo (D)
% A matriz dat.coeff_F dá os pesos de cada função objetiva

% Checagem de erros
dat = error_check_cst_3D_SPECIFY(dat); 

% Template dos structs
% Forma da planta
empty.type = dat.type;
empty.b = [];
empty.b1 = []; 
empty.c_r = []; 
empty.c_m = []; 
empty.c_t = [];
empty.mac = [];
empty.S = [];
empty.sweep = [];
empty.sweep1 = [];
empty.sweep2 = [];
% Aerofólios 
empty.v_ex_r = zeros(1,dat.BPn_r+2);
empty.v_in_r = zeros(1,dat.BPn_r+2);
empty.symm_r = [];
empty.v_ex_m = zeros(1,dat.BPn_m+2);
empty.v_in_m = zeros(1,dat.BPn_m+2);
empty.symm_m = [];
empty.v_ex_t = zeros(1,dat.BPn_t+2);
empty.v_in_t = zeros(1,dat.BPn_t+2);
empty.symm_t = [];
empty.tw_m = [];
empty.tw_t = [];
% Dados da malha
empty.NODE = [];
empty.PANEL = [];
empty.sec_N = [];
% Dados aerodinâmicos e pontuação
empty.aero = []; % Terá o mesmo formato que a matriz coeff_op
empty.score = 0;

% Inicializar os structs
pop = repmat(empty,dat.N,1);
chi = pop;

% Fazer um struct da asa original
original = empty;
original.type = dat.type;
original.b = dat.or_b;
original.b1 = dat.or_b1;
original.c_r = dat.or_c_r;
original.c_m = dat.or_c_m;
original.c_t = dat.or_c_t;
original.v_ex_r = dat.or_v_ex_r;
original.sweep = dat.or_sweep;
original.sweep1 = dat.or_sweep1;
original.sweep2 = dat.or_sweep2;
if dat.symm_override_r == 1
    original.v_in_r = dat.or_v_ex_r;
else
    original.v_in_r = dat.or_v_in_r;
end
original.v_ex_m = dat.or_v_ex_m;
if dat.symm_override_m == 1
    original.v_in_m = dat.or_v_ex_m;
else
    original.v_in_m = dat.or_v_in_m;
end
original.v_ex_t = dat.or_v_ex_t;
if dat.symm_override_t == 1
    original.v_in_t = dat.or_v_ex_t;
else
    original.v_in_t = dat.or_v_in_t;
end
original.tw_m = 'L';
original.tw_t = dat.or_tw_t;
% Pegar dados aerodinâmicos da asa original para fins de comparação
if dat.half == 1 % Considerar apenas metade da asa
    original = run_apame_mesher_cst_half(original,dat);
else % Considerar a asa completa
    original = run_apame_mesher_cst(original,dat);
end
original.aero = run_apame(original,dat);
if original.aero == 'n'
    i = input('A simulação do indivíduo original não convergiu. Continuar mesmo assim? (y/n) ','s');
    if i ~= 'y' && i ~= 'Y'
        error('Execução terminada')
    end
end


% Gerar população inicial
for i = 1:dat.N
    disp(['Indivíduo ' num2str(i)])
    
    % Gerar forma da planta
	% Envergadura
    pop(i).b = dat.or_b + rand*dat.b_ext(randi(2));
    if pop(i).b < dat.b_ext(3)
        pop(i).b = dat.b_ext(3);
    end
	% Corda da raiz
    pop(i).c_r = dat.or_c_r + rand*dat.c_r_ext(randi(2));
    if pop(i).c_r < dat.c_r_ext(3)
        pop(i).c_r = dat.c_r_ext(3);
    end
	% Corda da ponta
    pop(i).c_t = dat.or_c_t + rand*dat.c_t_ext(randi(2));
    if pop(i).c_t < dat.c_t_ext(3) % Se for menor que o mínimo permitido
        pop(i).c_t = dat.c_t_ext(3);
    elseif pop(i).c_t > pop(i).c_r % Se for maior que a corda da raiz
        pop(i).c_t = pop(i).c_r;
    end
	% Torção geométrica na ponta
    pop(i).tw_t = dat.or_tw_t + rand*dat.tw_t_ext(randi(2));
    
    % Gerar aerofólio da raiz
    check = 0;
    while check == 0
        
        if isequal(dat.or_v_ex_r,dat.or_v_in_r) || dat.symm_override_r == 1 % Perfil simétrico
            % Vetor com informações do extradorso
            pop(i).v_ex_r(1) = dat.or_v_ex_r(1) + rand*dat.le_R_ext1_r(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn_r
                pop(i).v_ex_r(a) = dat.or_v_ex_r(a) + rand*dat.A_ext1_r(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex_r(dat.BPn_r+1) = dat.or_v_ex_r(dat.BPn_r+1) + rand*dat.B_ext1_r(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex_r(dat.BPn_r+2) = dat.or_v_ex_r(dat.BPn_r+2); %randi(dat.delta_range)*0.1*rand;   % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in_r = pop(i).v_ex_r;
        
            % Checar o raio do bordo de ataque
            if pop(i).v_ex_r(1) < dat.le_R_ext1_r(3)
                pop(i).v_ex_r(1) = dat.le_R_ext1_r(3);
                pop(i).v_in_r(1) = pop(i).v_ex_r(1);
            end
            
            % Checar separação do bordo de fuga
            if pop(i).v_ex_r(dat.BPn_r+1) + pop(i).v_in_r(dat.BPn_r+1) < dat.B_ext2_r(1)
                pop(i).v_ex_r(dat.BPn_r+1) = dat.B_ext2_r(1)/2;
                pop(i).v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1)/2;
            end
            
            pop(i).symm_r = 1;
            
        else % perfil assimétrico
            % Vetor com informações do extradorso
            pop(i).v_ex_r(1) = dat.or_v_ex_r(1) + rand*dat.le_R_ext1_r(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn_r
                pop(i).v_ex_r(a) = dat.or_v_ex_r(a) + rand*dat.A_ext1_r(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex_r(dat.BPn_r+1) = dat.or_v_ex_r(dat.BPn_r+1) + rand*dat.B_ext1_r(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex_r(dat.BPn_r+2) = dat.or_v_ex_r(dat.BPn_r+2);  % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in_r(1) = dat.or_v_in_r(1) + rand*dat.le_R_ext2_r(randi(2));                     % Raio do bordo de ataque
            for a = 2:dat.BPn_r
                pop(i).v_in_r(a) = dat.or_v_in_r(a) + rand*dat.A_ext2_r(randi(2));        % Pesos intermediários
            end
            pop(i).v_in_r(dat.BPn_r+1) = (dat.B_ext2_r(1) - pop(i).v_ex_r(dat.BPn_r+1)) + rand*dat.B_ext2_r(2);     % Ângulo do bordo de fuga
            pop(i).v_in_r(dat.BPn_r+2) = dat.or_v_in_r(dat.BPn_r+2); 
            
            % Checar o raio do bordo de ataque
            if pop(i).v_ex_r(1) < dat.le_R_ext1_r(3)
                pop(i).v_ex_r(1) = dat.le_R_ext1_r(3);
            end
            if pop(i).v_in_r(1) < dat.le_R_ext2_r(3)
                pop(i).v_in_r(1) = dat.le_R_ext2_r(3);
            end
            
            % Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop(i).v_in_r(dat.BPn_r+1) < dat.B_ext2_r(1) - pop(i).v_ex_r(dat.BPn_r+1)
                pop(i).v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1) - pop(i).v_ex_r(dat.BPn_r+1);
            end
            
            % Checar os pesos (soma de pesos do intradorso deve ser menor ou
            % igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop(i).v_ex_r(2:dat.BPn_r));
            sum2 = sum(pop(i).v_in_r(2:dat.BPn_r));
            if sum2 > sum1,continue,end
            
            pop(i).symm_r = 0;
        end
        
        % Checagem de qualidade
        check = quality(run_cst_TCC2_3D(pop(i).v_ex_r,pop(i).v_in_r,dat,[dat.N1_r,dat.N2_r]),dat);
        
    end
    
    % Gerar aerofólio da ponta
    check = 0;
    while check == 0
        
        if isequal(dat.or_v_ex_t,dat.or_v_in_t) || dat.symm_override_t == 1 % Perfil simétrico
            % Vetor com informações do extradorso
            pop(i).v_ex_t(1) = dat.or_v_ex_t(1) + rand*dat.le_R_ext1_t(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn_t
                pop(i).v_ex_t(a) = dat.or_v_ex_t(a) + rand*dat.A_ext1_t(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex_t(dat.BPn_t+1) = dat.or_v_ex_t(dat.BPn_t+1) + rand*dat.B_ext1_t(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex_t(dat.BPn_t+2) = dat.or_v_ex_t(dat.BPn_t+2); %randi(dat.delta_tange)*0.1*rand;   % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in_t = pop(i).v_ex_t;
        
            % Checar o raio do bordo de ataque
            if pop(i).v_ex_t(1) < dat.le_R_ext1_t(3)
                pop(i).v_ex_t(1) = dat.le_R_ext1_t(3);
                pop(i).v_in_t(1) = pop(i).v_ex_t(1);
            end
            
            % Checar separação do bordo de fuga
            if pop(i).v_ex_t(dat.BPn_t+1) + pop(i).v_in_t(dat.BPn_t+1) < dat.B_ext2_t(1)
                pop(i).v_ex_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
                pop(i).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
            end
            
            pop(i).symm_t = 1;
            
        else % perfil assimétrico
            % Vetor com informações do extradorso
            pop(i).v_ex_t(1) = dat.or_v_ex_t(1) + rand*dat.le_R_ext1_t(randi(2));    % Raio do bordo de ataque
            for a = 2:dat.BPn_t
                pop(i).v_ex_t(a) = dat.or_v_ex_t(a) + rand*dat.A_ext1_t(randi(2));        % Pesos intermediários
            end
            pop(i).v_ex_t(dat.BPn_t+1) = dat.or_v_ex_t(dat.BPn_t+1) + rand*dat.B_ext1_t(randi(2));     % Ângulo do bordo de fuga
            pop(i).v_ex_t(dat.BPn_t+2) = dat.or_v_ex_t(dat.BPn_t+2);  % delta_z
            
            % Vetor com informações do intradorso
            pop(i).v_in_t(1) = dat.or_v_in_t(1) + rand*dat.le_R_ext2_t(randi(2));                     % Raio do bordo de ataque
            for a = 2:dat.BPn_t
                pop(i).v_in_t(a) = dat.or_v_in_t(a) + rand*dat.A_ext2_t(randi(2));        % Pesos intermediários
            end
            pop(i).v_in_t(dat.BPn_t+1) = (dat.B_ext2_t(1) - pop(i).v_ex_t(dat.BPn_t+1)) + rand*dat.B_ext2_t(2);     % Ângulo do bordo de fuga
            pop(i).v_in_t(dat.BPn_t+2) = dat.or_v_in_t(dat.BPn_t+2); 
            
            % Checar o raio do bordo de ataque
            if pop(i).v_ex_t(1) < dat.le_R_ext1_t(3)
                pop(i).v_ex_t(1) = dat.le_R_ext1_t(3);
            end
            if pop(i).v_in_t(1) < dat.le_R_ext2_t(3)
                pop(i).v_in_t(1) = dat.le_R_ext2_t(3);
            end
            
            % Checar separação do bordo de fuga (beta2>=L-beta1)
            if pop(i).v_in_t(dat.BPn_t+1) < dat.B_ext2_t(1) - pop(i).v_ex_t(dat.BPn_t+1)
                pop(i).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1) - pop(i).v_ex_t(dat.BPn_t+1);
            end
            
            % Checar os pesos (soma de pesos do intradorso deve ser menor ou
            % igual à soma de pesos do extradorso - pesos intermediários)
            sum1 = sum(pop(i).v_ex_t(2:dat.BPn_t));
            sum2 = sum(pop(i).v_in_t(2:dat.BPn_t));
            if sum2 > sum1,continue,end
            
            pop(i).symm_t = 0;
        end
        
        % Checagem de qualidade
        check = quality(run_cst_TCC2_3D(pop(i).v_ex_t,pop(i).v_in_t,dat,[dat.N1_t,dat.N2_t]),dat);
        
    end
    
    % Dados adicionais para asas bitrapezoidais
    if pop(i).type == 1
        % Mais dados da planta
        if dat.or_c_m == 'L'
            % Caso c_m seja definido como 'L', seu valor real será atribúido pela
            % função run_apame_mesher_cst
            pop(i).c_m = 'L';
        else    
            % Caso contrário, é realizado o processo abaixo
            pop(i).c_m = dat.or_c_m + rand*dat.c_m_ext(randi(2));
            if pop(i).c_m < dat.c_m_ext(3)
                pop(i).c_m = dat.c_m_ext(3);
            end
        end
        
        pop(i).b1 = dat.or_b1 + rand*dat.b1_ext(randi(2));
        if pop(i).b1 < dat.b1_ext(3)
            pop(i).b1 = dat.b1_ext(3);
        end
        
        % Aerofólio do meio
        pop(i).tw_m = 'L'; % Essa sempre será a configuração deste algoritmo, mas isso pode ser alterado (com as devidas alterações no resto do código)
        check = 0;
        while check == 0
            
            if isequal(dat.or_v_ex_m,dat.or_v_in_m) || dat.symm_override_m == 1 % Perfil simétrico
                % Vetor com informações do extradorso
                pop(i).v_ex_m(1) = dat.or_v_ex_m(1) + rand*dat.le_R_ext1_m(randi(2));    % Raio do bordo de ataque
                for a = 2:dat.BPn_m
                    pop(i).v_ex_m(a) = dat.or_v_ex_m(a) + rand*dat.A_ext1_m(randi(2));        % Pesos intermediários
                end
                pop(i).v_ex_m(dat.BPn_m+1) = dat.or_v_ex_m(dat.BPn_m+1) + rand*dat.B_ext1_m(randi(2));     % Ângulo do bordo de fuga
                pop(i).v_ex_m(dat.BPn_m+2) = dat.or_v_ex_m(dat.BPn_m+2); %randi(dat.delta_mange)*0.1*rand;   % delta_z
                
                % Vetor com informações do intradorso
                pop(i).v_in_m = pop(i).v_ex_m;
            
                % Checar o raio do bordo de ataque
                if pop(i).v_ex_m(1) < dat.le_R_ext1_m(3)
                    pop(i).v_ex_m(1) = dat.le_R_ext1_m(3);
                    pop(i).v_in_m(1) = pop(i).v_ex_m(1);
                end
                
                % Checar separação do bordo de fuga
                if pop(i).v_ex_m(dat.BPn_m+1) + pop(i).v_in_m(dat.BPn_m+1) < dat.B_ext2_m(1)
                    pop(i).v_ex_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
                    pop(i).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
                end
                
                pop(i).symm_m = 1;
                
            else % perfil assimétrico
                % Vetor com informações do extradorso
                pop(i).v_ex_m(1) = dat.or_v_ex_m(1) + rand*dat.le_R_ext1_m(randi(2));    % Raio do bordo de ataque
                for a = 2:dat.BPn_m
                    pop(i).v_ex_m(a) = dat.or_v_ex_m(a) + rand*dat.A_ext1_m(randi(2));        % Pesos intermediários
                end
                pop(i).v_ex_m(dat.BPn_m+1) = dat.or_v_ex_m(dat.BPn_m+1) + rand*dat.B_ext1_m(randi(2));     % Ângulo do bordo de fuga
                pop(i).v_ex_m(dat.BPn_m+2) = dat.or_v_ex_m(dat.BPn_m+2);  % delta_z
                
                % Vetor com informações do intradorso
                pop(i).v_in_m(1) = dat.or_v_in_m(1) + rand*dat.le_R_ext2_m(randi(2));                     % Raio do bordo de ataque
                for a = 2:dat.BPn_m
                    pop(i).v_in_m(a) = dat.or_v_in_m(a) + rand*dat.A_ext2_m(randi(2));        % Pesos intermediários
                end
                pop(i).v_in_m(dat.BPn_m+1) = (dat.B_ext2_m(1) - pop(i).v_ex_m(dat.BPn_m+1)) + rand*dat.B_ext2_m(2);     % Ângulo do bordo de fuga
                pop(i).v_in_m(dat.BPn_m+2) = dat.or_v_in_m(dat.BPn_m+2); 
                
                % Checar o raio do bordo de ataque
                if pop(i).v_ex_m(1) < dat.le_R_ext1_m(3)
                    pop(i).v_ex_m(1) = dat.le_R_ext1_m(3);
                end
                if pop(i).v_in_m(1) < dat.le_R_ext2_m(3)
                    pop(i).v_in_m(1) = dat.le_R_ext2_m(3);
                end
                
                % Checar separação do bordo de fuga (beta2>=L-beta1)
                if pop(i).v_in_m(dat.BPn_m+1) < dat.B_ext2_m(1) - pop(i).v_ex_m(dat.BPn_m+1)
                    pop(i).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1) - pop(i).v_ex_m(dat.BPn_m+1);
                end
                
                % Checar os pesos (soma de pesos do intradorso deve ser menor ou
                % igual à soma de pesos do extradorso - pesos intermediários)
                sum1 = sum(pop(i).v_ex_m(2:dat.BPn_m));
                sum2 = sum(pop(i).v_in_m(2:dat.BPn_m));
                if sum2 > sum1,continue,end
                
                pop(i).symm_m = 0;
            end
            
            % Checagem de qualidade
            check = quality(run_cst_TCC2_3D(pop(i).v_ex_m,pop(i).v_in_m,dat,[dat.N1_m,dat.N2_m]),dat);
            
        end
       
        % Aplicar checagens de geometria da planta
        if pop(i).c_m > pop(i).c_r && dat.or_c_m ~= 'L'
            pop(i).c_m = pop(i).c_r;
        end
        if pop(i).c_t > pop(i).c_m && dat.or_c_m ~= 'L'
            pop(i).c_t = pop(i).c_m;
        end
        if pop(i).b1 > pop(i).b - 2*dat.b1_ext(4)
            pop(i).b1 = pop(i).b - 2*dat.b1_ext(4);
        end
    
        % Enflechamento da primeira seção
        if dat.or_sweep1 == 'Z'
            pop(i).sweep1 = 'Z';
        else
            pop(i).sweep1 = dat.or_sweep1 + rand*dat.sweep1_ext(randi(2));
        end
        
        % Enflechamento da segunda seção
        if  dat.or_sweep2 == 'Z'
            pop(i).sweep2 = 'Z';
        else
            pop(i).sweep2 = dat.or_sweep2 + rand*dat.sweep2_ext(randi(2));
        end
        
        
        % Debugging: ver se a opção 'L' da corda do meio está sendo cumprida
        if dat.or_c_m == 'L' && pop(i).c_m ~= 'L'
            error("Erro da opção 'L' da corda do meio")
        end
        
        
    else % Estabelecer enflechamento da asa trapezoidal simples 
        if dat.or_sweep == 'Z'
            pop(i).sweep = 'Z';
        else
            pop(i).sweep = dat.or_sweep + rand*dat.sweep_ext(randi(2));
        end
        
    end
    
end

% Gerar struct que guarda o melhor perfil de cada geração
archive = repmat(empty,dat.iter,1);

%% Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for loop = 1:dat.iter
    
    % Simular as asas e obter dados
    select = ones(1,dat.N);
    disp('<< Simulação das asas >>')
    for i = 1:dat.N
        disp(['Indivíduo ' num2str(i)])
        if dat.half == 1 % Considerar apenas metade da asa
            pop(i) = run_apame_mesher_cst_half(pop(i),dat); % Obter malha
        else % Considerar a asa completa
            pop(i) = run_apame_mesher_cst(pop(i),dat); % Obter malha
        end
        pop(i).aero = run_apame(pop(i),dat); % Simulação
        
        % Marcar indivíduos que não tenham convergido na simulação
        if pop(i).aero == 'n' || isempty(pop(i).aero)
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
        error('Nenhuma asa convergiu nas simulações')
	end
    
    % Atribuir pontuações (fitnesses)
    pop = fitness_cst_3D(pop,dat,select2);
    
    % Isto serve pra pôr todas as pontuações em um vetor
    weights = [pop.score];
    
    % Guardar a melhor asa de cada iteração
    [~,pos] = max(weights);
    archive(loop) = pop(pos);
%    for i = 1:4 % [apagar isto depois]
%        % As médias contabilizam apenas os indivíduos que convergiram nas simulações
%        temp = make_vector(pop,i,select2);
%        dat.aero_M(loop,i) = mean(temp(find(temp~=0)));
%    end
    
    % Mostrar a melhor asa e comparar com a original
%    figure(1),clf,hold on,figure(2),clf,hold on
%    text = ['Iteração ' num2str(loop)];
    figure(1),clf,hold on
    plot_planform(original,dat.half,'k')
    plot_planform(pop(pos),dat.half,'r--')
    title(['Iteração ' num2str(loop)])
    legend('Original','Otimizado')
%    run_apame_mesher_cst(pop(pos),dat,2,2,1,2,text);
%    run_apame_mesher_cst(original,dat,2,2,1,2,text);
    
    
    
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
		
        % Isto seleciona dois pais por meio de uma seleção via roleta
        % (indivíduos com pesos maiores têm  mais chance de serem selecionados)
        par = [0,0]; % Vetor que indica a numeração dos pais escolhidos
        par(1) = selection_crossover(weights);
        par(2) = selection_crossover(weights);
        
        % Crossover, gerando dois filhos de cada dois pais
        % Genes a serem trocados
        % Envergadura b
        % Envergadura b1 (asas bitrapezoidais apenas)
        % Corda da raiz
        % Corda do meio (asas bitrapezoidais apenas)
        % Corda da ponta
        % - aerofólio da raiz
        % - aerofólio do meio (asas bitrapezoidais apenas)
        % - aerofólio da ponta
        % Torção geométrica na ponta
		% Enflechamento total (asas trapezoidais simples)
        % Enflechamento da primeira seção (asas bitrapezoidais)
        % Enflechamento da segunda seção (asas bitrapezoidais)
        
%        
%        if pop(par(1)).type ~= pop(par(2)).type % Caso as asas sejam de tipos diferentes
%            op = randi([1,2]);
%        else
%            op = 1; % Caso as asas sejam de tipos iguais
%        end
        
        op = 1;
        
        if op == 1 % Fazer crossover normalmente (trocar genes gerais)
            % Escrever um vetor de opções:
            % [b,b1,c_r,c_m,c_t,airfoil_r,airfoil_m,airfoil_t,tw_t]
            cross_op_v = randi([0,1],1,12);
            
            chi(c) = pop(par(1));
            chi(c+1) = pop(par(2));
            
            % Envergadura b
            if cross_op_v(1) == 1
                chi(c).b = pop(par(2)).b;                
                chi(c+1).b = pop(par(1)).b;
            end 
        
            % Envergadura da primeira seção b1 (asas bitrapezoidais apenas) (deve cumprir o requisito de separação mínima)
            if cross_op_v(2) == 1 && pop(par(1)).type == 1 && pop(par(2)).type == 1 %&& pop(par(1)).b1 <= pop(par(2)).b - dat.b1_ext_in(3)*2 && pop(par(2)).b1 <= pop(par(1)).b - dat.b1_ext_in(3)*2
                chi(c).b1 = pop(par(2)).b1;
                chi(c+1).b1 = pop(par(1)).b1;
            end
            
            % Corda da raiz c_r (deve cumprir o requisito c_r >= c_t)
%            if cross_op_v(3) == 1 && pop(par(1)).c_r >= pop(par(2)).c_t && pop(par(2)).c_r >= pop(par(1)).c_t
            if cross_op_v(3) == 1 && chi(c).c_r >= chi(c+1).c_t && chi(c+1).c_r >= chi(c).c_t
                chi(c).c_r = pop(par(2)).c_r;
                chi(c+1).c_r = pop(par(1)).c_r;
            end
            
            % Corda da ponta c_t (deve cumprir o requisito c_r >= c_t)
%            if cross_op_v(4) == 1 && pop(par(1)).c_r >= pop(par(2)).c_t && pop(par(2)).c_r >= pop(par(1)).c_t 
            if cross_op_v(4) == 1 && chi(c).c_r >= chi(c+1).c_t && chi(c+1).c_r >= chi(c).c_t
                chi(c).c_t = pop(par(2)).c_t;
                chi(c+1).c_t = pop(par(1)).c_t;
            end
            
            % Corda do meio c_m (asas bitrapezoidais apenas) (deve cumprir o requisito c_r >= c_m >= c_t) (ignorar caso a opção 'L' seja aplicada a c_m)
            if cross_op_v(5) == 1 && pop(par(1)).type == 1 && pop(par(2)).type == 1 && ~dat.or_c_m == 'L' %&& pop(par(1)).c_r >= pop(par(2)).c_m && pop(par(1)).c_m >= pop(par(2)).c_t && pop(par(2)).c_r >= pop(par(1)).c_m && pop(par(2)).c_m >= pop(par(1)).c_t 
                chi(c).c_m = pop(par(2)).c_m;
                chi(c+1).c_m = pop(par(1)).c_m;
            end
            
            % Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi(c).type == 1
                
                % Se a envergadura da primeira seção for maior do que permitido
                % pelo requisito de separação, atribuir o máximo valor que 
                % cumpre o requisito
                if chi(c).b1 > chi(c).b-dat.b1_ext(4)*2
                    chi(c).b1 = chi(c).b-dat.b1_ext(4)*2;
                end
                
                % Se a corda do meio for maior do que a corda da raiz, atribuir
                % o valor da raiz ao meio
                if chi(c).c_m > chi(c).c_r && dat.or_c_m ~= 'L'
                    chi(c).c_m = chi(c).c_r;
                end
                
                % Se a corda do meio for menor que a corda da ponta, atribuir
                % o valor da ponta ao meio
                if chi(c).c_m < chi(c).c_t && dat.or_c_m ~= 'L'
                    chi(c).c_m = chi(c).c_t;
                end
                
            end
            if chi(c+1).type == 1
                % Se a envergadura da primeira seção for maior do que permitido
                % pelo requisito de separação, atribuir o máximo valor que 
                % cumpre o requisito
                if chi(c+1).b1 > chi(c+1).b-dat.b1_ext(4)*2
                    chi(c+1).b1 = chi(c+1).b-dat.b1_ext(4)*2;
                end
                
                % Se a corda do meio for maior do que a corda da raiz, atribuir
                % o valor da raiz ao meio
                if chi(c+1).c_m > chi(c+1).c_r && dat.or_c_m ~= 'L'
                    chi(c+1).c_m = chi(c+1).c_r;
                end
                
                % Se a corda do meio for menor que a corda da ponta, atribuir
                % o valor da ponta ao meio
                if chi(c+1).c_m < chi(c+1).c_t && dat.or_c_m ~= 'L'
                    chi(c+1).c_m = chi(c+1).c_t;
                end
            end
            
            % Aerofólio da raiz
            if cross_op_v(6) == 1 
            
                af1_ex = chi(c).v_ex_r;
                af1_in = chi(c).v_in_r;
                af2_ex = chi(c+1).v_ex_r;
                af2_in = chi(c+1).v_in_r;
                
                check = 0;
                while check == 0;
                
                    m = randi([1,2]);
                    switch m
                        case 1 % Trocar os perfis completamente
                            chi(c).v_ex_r = af2_ex;
                            chi(c).v_in_r = af2_in;
                            chi(c+1).v_ex_r = af1_ex;
                            chi(c+1).v_in_r = af1_in;
                        
                        case 2 % Fazer o crossover das características
                            
                            % v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                            if chi(c).symm_r == 0 && chi(c+1).symm_r == 0 % Se ambos forem assimétricos
                                n = randi([1,4]);
                                switch n
                                    case 1 % Trocar os extradorsos e intradorsos inteiros
                                        temp1_ex = af1_ex;
                                        temp1_in = af2_in;
                                        temp2_ex = af2_ex;
                                        temp2_in = af1_in;
                                        
                                    case 2 % Trocar o raio do bordo de ataque
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = [af1_in(1),af2_in(2:end)];
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = [af2_in(1),af1_in(2:end)];
                                        
                                    case 3 % Trocar os pesos intermediários
                                        if dat.BPn_r == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_r),af2_ex(dat.BPn_r+1:end)];
                                            temp1_in = [af2_in(1),af1_in(2:dat.BPn_r),af2_in(dat.BPn_r+1:end)];
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_r),af1_ex(dat.BPn_r+1:end)];
                                            temp2_in = [af1_in(1),af2_in(2:dat.BPn_r),af1_in(dat.BPn_r+1:end)];
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_r);
                                            num2 = randi(2:dat.BPn_r);

                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_r)];
                                            temp1_2 = [af1_in(2:num2),af2_in(num2+1:dat.BPn_r)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_r)];
                                            temp2_2 = [af2_in(2:num2),af1_in(num2+1:dat.BPn_r)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_r+1:end)];
                                            temp1_in = [af1_in(1),temp2_2,af1_in(dat.BPn_r+1:end)];
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_r+1:end)];
                                            temp2_in = [af2_in(1),temp1_2,af2_in(dat.BPn_r+1:end)]; 
                                            
                                        end
                                        
                                    case 4 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_r),af1_ex(dat.BPn_r+1),af2_ex(dat.BPn_r+2)];
                                        temp1_in = [af2_in(1:dat.BPn_r),af1_in(dat.BPn_r+1),af2_in(dat.BPn_r+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_r),af2_ex(dat.BPn_r+1),af1_ex(dat.BPn_r+2)];
                                        temp2_in = [af1_in(1:dat.BPn_r),af2_in(dat.BPn_r+1),af1_in(dat.BPn_r+2)];
                                            
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_r = temp1_ex;
                                    chi(c).v_in_r = temp1_in;
                                    chi(c+1).v_ex_r = temp2_ex;
                                    chi(c+1).v_in_r = temp2_in;
                                else
                                    chi(c).v_ex_r = temp2_ex;
                                    chi(c).v_in_r = temp2_in;
                                    chi(c+1).v_ex_r = temp1_ex;
                                    chi(c+1).v_in_r = temp1_in;
                                end
                                
                                if n == 1
                                        %             chi(c) = pop(par(1)); af1  [(referência)]
                                    %            chi(c+1) = pop(par(2)); af2 
                                    
                                    
                                    % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                    % requisito de separação, alterar o ângulo do bordo de fuga
                                    % do intradorso
                                    if chi(c).v_in_r(dat.BPn_r+1) < (dat.B_ext2_r(1)-chi(c).v_ex_r(dat.BPn_r+1))
                                        chi(c).v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1)-chi(c).v_ex_r(dat.BPn_r+1);
                                    end
                                    if chi(c+1).v_in_r(dat.BPn_r+1) < (dat.B_ext2_r(1)-chi(c+1).v_ex_r(dat.BPn_r+1))
                                        chi(c+1).v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1)-chi(c+1).v_ex_r(dat.BPn_r+1);
                                    end
                                    % Decisão de alterar o intradorso em base de uma nota na página
                                    % 57(87) do Raymer (2018)
                                end
                                
                                % Checar os pesos
                                sum1 = sum(chi(c).v_ex_r(2:dat.BPn_r));
                                sum2 = sum(chi(c).v_in_r(2:dat.BPn_r));
                                if sum2 > sum1,continue,end
                                sum1 = sum(chi(c+1).v_ex_r(2:dat.BPn_r));
                                sum2 = sum(chi(c+1).v_in_r(2:dat.BPn_r));
                                if sum2 > sum1,continue,end
                                
                                % Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                                % for usado como variável)
                                %chi(c).v_in(dat.BPn) = -chi(c).v_ex(dat.BPn)*chi(c).v_ex(dat.BPn-1)/chi(c).v_in(dat.BPn-1);
                                %chi(c+1).v_in(dat.BPn) = -chi(c+1).v_ex(dat.BPn)*chi(c+1).v_ex(dat.BPn-1)/chi(c+1).v_in(dat.BPn-1);
                            
                                chi(c).symm_r = 0;
                                chi(c+1).symm_r = 0;
                        
                            elseif chi(c).symm_r == 1 && chi(c+1).symm_r == 1 % Se ambos forem simétricos
                                n = randi([1,3]);
                                switch n                        
                                    case 1 % Trocar o raio do bordo de ataque
                                        
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = temp1_ex;
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = temp2_ex;
                                        
                                    case 2 % Trocar os pesos intermediários
                                        if dat.BPn_r == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_r),af2_ex(dat.BPn_r+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_r),af1_ex(dat.BPn_r+1:end)];
                                            temp2_in = temp2_ex;
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_r);
                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_r)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_r)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_r+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_r+1:end)];
                                            temp2_in = temp2_ex;

                                        end
                                        
                                    case 3 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_r),af1_ex(dat.BPn_r+1),af2_ex(dat.BPn_r+2)];
                                        temp1_in = [af2_in(1:dat.BPn_r),af1_in(dat.BPn_r+1),af2_in(dat.BPn_r+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_r),af2_ex(dat.BPn_r+1),af1_ex(dat.BPn_r+2)];
                                        temp2_in = [af1_in(1:dat.BPn_r),af2_in(dat.BPn_r+1),af1_in(dat.BPn_r+2)];
                                        
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_r = temp1_ex;
                                    chi(c).v_in_r = temp1_in;
                                    chi(c+1).v_ex_r = temp2_ex;
                                    chi(c+1).v_in_r = temp2_in;
                                else
                                    chi(c).v_ex_r = temp2_ex;
                                    chi(c).v_in_r = temp2_in;
                                    chi(c+1).v_ex_r = temp1_ex;
                                    chi(c+1).v_in_r = temp1_in;
                                end
                                
                                chi(c).symm_r = 1;
                                chi(c+1).symm_r = 1;
                            
                            else % Se um for simétrico e o outro for assimétrico 
                                % Transformar o simétrico em um assimétrico e vice-versa
                                if chi(c).symm_r == 0 
                                    temp1_ex = af1_ex;
                                    temp1_in = af1_ex;
                                    temp2_ex = af2_ex;
                                    temp2_in = af1_in;
                                    
                                else
                                    temp1_ex = af2_ex;
                                    temp1_in = af2_ex;
                                    temp2_ex = af1_ex;
                                    temp2_in = af2_in;
                                    
                                end

                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_r = temp1_ex;
                                    chi(c).v_in_r = temp1_in;
                                    chi(c+1).v_ex_r = temp2_ex;
                                    chi(c+1).v_in_r = temp2_in;
                                else
                                    chi(c).v_ex_r = temp2_ex;
                                    chi(c).v_in_r = temp2_in;
                                    chi(c+1).v_ex_r = temp1_ex;
                                    chi(c+1).v_in_r = temp1_in;
                                end
                                
                                % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                % requisito de separação, alterar o ângulo do bordo de fuga
                                % do intradorso
                                if chi(c).v_in(dat.BPn_r+1) < (dat.B_ext2_r(1)-chi(c).v_ex(dat.BPn_r+1))
                %                    chi(c).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1);
                                    chi(c).v_ex(dat.BPn_r+1) = dat.B_ext2(1)/2;
                                    chi(c).v_in(dat.BPn_r+1) = dat.B_ext2(1)/2;
                                end
                                if chi(c+1).v_in(dat.BPn_r+1) < (dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn_r+1))
                                    chi(c+1).v_in(dat.BPn_r+1) = dat.B_ext2(1)-chi(c+1).v_ex(dat.BPn_r+1);
                                end
                                % Decisão de alterar o intradorso em base de uma nota na página
                                % 57(87) do Raymer (2018)
                                
                                % Checar os pesos
                                sum1 = sum(chi(c+1).v_ex_r(2:dat.BPn_r));
                                sum2 = sum(chi(c+1).v_in_r(2:dat.BPn_r));
                                if sum2 > sum1,continue,end
                                
                                chi(c).symm_r = 1;
                                chi(c+1).symm_r = 0;
                                
                            end
                    end 
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi(c).v_ex_r,chi(c).v_in_r,dat,[dat.N1_r,dat.N2_r]),dat);
                    if check == 0,continue,end
                    check = quality(run_cst_TCC2_3D(chi(c+1).v_ex_r,chi(c+1).v_in_r,dat,[dat.N1_r,dat.N2_r]),dat); 
                
                end
            end
            
            % Aerofólio do meio (asas bitrapezoidais apenas) (ignorar caso a opção 'L' seja aplicada ao aerofólio do meio)
            if cross_op_v(7) == 1 && pop(par(1)).type == 1 && pop(par(2)).type == 1 && ~ismember('L',dat.le_R_ext1_m)
                af1_ex = chi(c).v_ex_m;
                af1_in = chi(c).v_in_m;
                af2_ex = chi(c+1).v_ex_m;
                af2_in = chi(c+1).v_in_m;
                
                check = 0;
                while check == 0;
                
                    m = randi([1,2]);
                    switch m
                        case 1 % Trocar os perfis completamente
                            chi(c).v_ex_m = af2_ex;
                            chi(c).v_in_m = af2_in;
                            chi(c+1).v_ex_m = af1_ex;
                            chi(c+1).v_in_m = af1_in;
                        
                        case 2 % Fazer o crossover das características
                            
                            % v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                            if chi(c).symm_m == 0 && chi(c+1).symm_m == 0 % Se ambos forem assimétricos
                                n = randi([1,4]);
                                switch n
                                    case 1 % Trocar os extradorsos e intradorsos inteiros
                                        temp1_ex = af1_ex;
                                        temp1_in = af2_in;
                                        temp2_ex = af2_ex;
                                        temp2_in = af1_in;
                                        
                                    case 2 % Trocar o raio do bordo de ataque
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = [af1_in(1),af2_in(2:end)];
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = [af2_in(1),af1_in(2:end)];
                                        
                                    case 3 % Trocar os pesos intermediários
                                        if dat.BPn_m == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_m),af2_ex(dat.BPn_m+1:end)];
                                            temp1_in = [af2_in(1),af1_in(2:dat.BPn_m),af2_in(dat.BPn_m+1:end)];
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_m),af1_ex(dat.BPn_m+1:end)];
                                            temp2_in = [af1_in(1),af2_in(2:dat.BPn_m),af1_in(dat.BPn_m+1:end)];
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_m);
                                            num2 = randi(2:dat.BPn_m);

                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_m)];
                                            temp1_2 = [af1_in(2:num2),af2_in(num2+1:dat.BPn_m)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_m)];
                                            temp2_2 = [af2_in(2:num2),af1_in(num2+1:dat.BPn_m)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_m+1:end)];
                                            temp1_in = [af1_in(1),temp2_2,af1_in(dat.BPn_m+1:end)];
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_m+1:end)];
                                            temp2_in = [af2_in(1),temp1_2,af2_in(dat.BPn_m+1:end)]; 
                                            
                                        end
                                        
                                    case 4 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_m),af1_ex(dat.BPn_m+1),af2_ex(dat.BPn_m+2)];
                                        temp1_in = [af2_in(1:dat.BPn_m),af1_in(dat.BPn_m+1),af2_in(dat.BPn_m+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_m),af2_ex(dat.BPn_m+1),af1_ex(dat.BPn_m+2)];
                                        temp2_in = [af1_in(1:dat.BPn_m),af2_in(dat.BPn_m+1),af1_in(dat.BPn_m+2)];
                                            
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_m = temp1_ex;
                                    chi(c).v_in_m = temp1_in;
                                    chi(c+1).v_ex_m = temp2_ex;
                                    chi(c+1).v_in_m = temp2_in;
                                else
                                    chi(c).v_ex_m = temp2_ex;
                                    chi(c).v_in_m = temp2_in;
                                    chi(c+1).v_ex_m = temp1_ex;
                                    chi(c+1).v_in_m = temp1_in;
                                end
                                
                                if n == 1
                                        %             chi(c) = pop(par(1)); af1  [(referência)]
                                    %            chi(c+1) = pop(par(2)); af2 
                                    
                                    
                                    % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                    % requisito de separação, alterar o ângulo do bordo de fuga
                                    % do intradorso
                                    if chi(c).v_in_m(dat.BPn_m+1) < (dat.B_ext2_m(1)-chi(c).v_ex_m(dat.BPn_m+1))
                                        chi(c).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)-chi(c).v_ex_m(dat.BPn_m+1);
                                    end
                                    if chi(c+1).v_in_m(dat.BPn_m+1) < (dat.B_ext2_m(1)-chi(c+1).v_ex_m(dat.BPn_m+1))
                                        chi(c+1).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)-chi(c+1).v_ex_m(dat.BPn_m+1);
                                    end
                                    % Decisão de alterar o intradorso em base de uma nota na página
                                    % 57(87) do Raymer (2018)
                                end
                                
                                % Checar os pesos
                                sum1 = sum(chi(c).v_ex_m(2:dat.BPn_m));
                                sum2 = sum(chi(c).v_in_m(2:dat.BPn_m));
                                if sum2 > sum1,continue,end
                                sum1 = sum(chi(c+1).v_ex_m(2:dat.BPn_m));
                                sum2 = sum(chi(c+1).v_in_m(2:dat.BPn_m));
                                if sum2 > sum1,continue,end
                                
                                % Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                                % for usado como variável)
                                %chi(c).v_in(dat.BPn) = -chi(c).v_ex(dat.BPn)*chi(c).v_ex(dat.BPn-1)/chi(c).v_in(dat.BPn-1);
                                %chi(c+1).v_in(dat.BPn) = -chi(c+1).v_ex(dat.BPn)*chi(c+1).v_ex(dat.BPn-1)/chi(c+1).v_in(dat.BPn-1);
                            
                                chi(c).symm_m = 0;
                                chi(c+1).symm_m = 0;
                        
                            elseif chi(c).symm_m == 1 && chi(c+1).symm_m == 1 % Se ambos forem simétricos
                                n = randi([1,3]);
                                switch n                        
                                    case 1 % Trocar o raio do bordo de ataque
                                        
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = temp1_ex;
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = temp2_ex;
                                        
                                    case 2 % Trocar os pesos intermediários
                                        if dat.BPn_m == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_m),af2_ex(dat.BPn_m+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_m),af1_ex(dat.BPn_m+1:end)];
                                            temp2_in = temp2_ex;
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_m);
                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_m)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_m)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_m+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_m+1:end)];
                                            temp2_in = temp2_ex;

                                        end
                                        
                                    case 3 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_m),af1_ex(dat.BPn_m+1),af2_ex(dat.BPn_m+2)];
                                        temp1_in = [af2_in(1:dat.BPn_m),af1_in(dat.BPn_m+1),af2_in(dat.BPn_m+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_m),af2_ex(dat.BPn_m+1),af1_ex(dat.BPn_m+2)];
                                        temp2_in = [af1_in(1:dat.BPn_m),af2_in(dat.BPn_m+1),af1_in(dat.BPn_m+2)];
                                        
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_m = temp1_ex;
                                    chi(c).v_in_m = temp1_in;
                                    chi(c+1).v_ex_m = temp2_ex;
                                    chi(c+1).v_in_m = temp2_in;
                                else
                                    chi(c).v_ex_m = temp2_ex;
                                    chi(c).v_in_m = temp2_in;
                                    chi(c+1).v_ex_m = temp1_ex;
                                    chi(c+1).v_in_m = temp1_in;
                                end
                                
                                chi(c).symm_m = 1;
                                chi(c+1).symm_m = 1;
                            
                            else % Se um for simétrico e o outro for assimétrico 
                                % Transformar o simétrico em um assimétrico e vice-versa
                                if chi(c).symm_m == 0 
                                    temp1_ex = af1_ex;
                                    temp1_in = af1_ex;
                                    temp2_ex = af2_ex;
                                    temp2_in = af1_in;
                                    
                                else
                                    temp1_ex = af2_ex;
                                    temp1_in = af2_ex;
                                    temp2_ex = af1_ex;
                                    temp2_in = af2_in;
                                    
                                end

                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_m = temp1_ex;
                                    chi(c).v_in_m = temp1_in;
                                    chi(c+1).v_ex_m = temp2_ex;
                                    chi(c+1).v_in_m = temp2_in;
                                else
                                    chi(c).v_ex_m = temp2_ex;
                                    chi(c).v_in_m = temp2_in;
                                    chi(c+1).v_ex_m = temp1_ex;
                                    chi(c+1).v_in_m = temp1_in;
                                end
                                
                                % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                % requisito de separação, alterar o ângulo do bordo de fuga
                                % do intradorso
                                if chi(c).v_in_m(dat.BPn_m+1) < (dat.B_ext2_m(1)-chi(c).v_ex_m(dat.BPn_m+1))
                %                    chi(c).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1);
                                    chi(c).v_ex_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
                                    chi(c).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
                                end
                                if chi(c+1).v_in_m(dat.BPn_m+1) < (dat.B_ext2_m(1)-chi(c+1).v_ex_m(dat.BPn_m+1))
                                    chi(c+1).v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)-chi(c+1).v_ex_m(dat.BPn_m+1);
                                end
                                % Decisão de alterar o intradorso em base de uma nota na página
                                % 57(87) do Raymer (2018)
                                
                                % Checar os pesos
                                sum1 = sum(chi(c+1).v_ex_m(2:dat.BPn_m));
                                sum2 = sum(chi(c+1).v_in_m(2:dat.BPn_m));
                                if sum2 > sum1,continue,end
                                
                                chi(c).symm_m = 1;
                                chi(c+1).symm_m = 0;
                                
                            end
                    end 
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi(c).v_ex_m,chi(c).v_in_m,dat,[dat.N1_m,dat.N2_m]),dat);
                    if check == 0,continue,end
                    check = quality(run_cst_TCC2_3D(chi(c+1).v_ex_m,chi(c+1).v_in_m,dat,[dat.N1_m,dat.N2_m]),dat); 
                
                end
            end
            
            % Aerofólio da ponta
            if cross_op_v(8) == 1 
                af1_ex = chi(c).v_ex_t;
                af1_in = chi(c).v_in_t;
                af2_ex = chi(c+1).v_ex_t;
                af2_in = chi(c+1).v_in_t;
                
                check = 0;
                while check == 0;
                
                    m = randi([1,2]);
                    switch m
                        case 1 % Trocar os perfis completamente
                            chi(c).v_ex_t = af2_ex;
                            chi(c).v_in_t = af2_in;
                            chi(c+1).v_ex_t = af1_ex;
                            chi(c+1).v_in_t = af1_in;
                        
                        case 2 % Fazer o crossover das características
                            
                            % v = [ RLe A1 A2 A3 ... A(N) beta Dz ]
                            if chi(c).symm_t == 0 && chi(c+1).symm_t == 0 % Se ambos forem assimétricos
                                n = randi([1,4]);
                                switch n
                                    case 1 % Trocar os extradorsos e intradorsos inteiros
                                        temp1_ex = af1_ex;
                                        temp1_in = af2_in;
                                        temp2_ex = af2_ex;
                                        temp2_in = af1_in;
                                        
                                    case 2 % Trocar o raio do bordo de ataque
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = [af1_in(1),af2_in(2:end)];
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = [af2_in(1),af1_in(2:end)];
                                        
                                    case 3 % Trocar os pesos intermediários
                                        if dat.BPn_t == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_t),af2_ex(dat.BPn_t+1:end)];
                                            temp1_in = [af2_in(1),af1_in(2:dat.BPn_t),af2_in(dat.BPn_t+1:end)];
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_t),af1_ex(dat.BPn_t+1:end)];
                                            temp2_in = [af1_in(1),af2_in(2:dat.BPn_t),af1_in(dat.BPn_t+1:end)];
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_t);
                                            num2 = randi(2:dat.BPn_t);

                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_t)];
                                            temp1_2 = [af1_in(2:num2),af2_in(num2+1:dat.BPn_t)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_t)];
                                            temp2_2 = [af2_in(2:num2),af1_in(num2+1:dat.BPn_t)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_t+1:end)];
                                            temp1_in = [af1_in(1),temp2_2,af1_in(dat.BPn_t+1:end)];
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_t+1:end)];
                                            temp2_in = [af2_in(1),temp1_2,af2_in(dat.BPn_t+1:end)]; 
                                            
                                        end
                                        
                                    case 4 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_t),af1_ex(dat.BPn_t+1),af2_ex(dat.BPn_t+2)];
                                        temp1_in = [af2_in(1:dat.BPn_t),af1_in(dat.BPn_t+1),af2_in(dat.BPn_t+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_t),af2_ex(dat.BPn_t+1),af1_ex(dat.BPn_t+2)];
                                        temp2_in = [af1_in(1:dat.BPn_t),af2_in(dat.BPn_t+1),af1_in(dat.BPn_t+2)];
                                            
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_t = temp1_ex;
                                    chi(c).v_in_t = temp1_in;
                                    chi(c+1).v_ex_t = temp2_ex;
                                    chi(c+1).v_in_t = temp2_in;
                                else
                                    chi(c).v_ex_t = temp2_ex;
                                    chi(c).v_in_t = temp2_in;
                                    chi(c+1).v_ex_t = temp1_ex;
                                    chi(c+1).v_in_t = temp1_in;
                                end
                                
                                if n == 1
                                        %             chi(c) = pop(par(1)); af1  [(referência)]
                                    %            chi(c+1) = pop(par(2)); af2 
                                    
                                    
                                    % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                    % requisito de separação, alterar o ângulo do bordo de fuga
                                    % do intradorso
                                    if chi(c).v_in_t(dat.BPn_t+1) < (dat.B_ext2_t(1)-chi(c).v_ex_t(dat.BPn_t+1))
                                        chi(c).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)-chi(c).v_ex_t(dat.BPn_t+1);
                                    end
                                    if chi(c+1).v_in_t(dat.BPn_t+1) < (dat.B_ext2_t(1)-chi(c+1).v_ex_t(dat.BPn_t+1))
                                        chi(c+1).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)-chi(c+1).v_ex_t(dat.BPn_t+1);
                                    end
                                    % Decisão de alterar o intradorso em base de uma nota na página
                                    % 57(87) do Raymer (2018)
                                end
                                
                                % Checar os pesos
                                sum1 = sum(chi(c).v_ex_t(2:dat.BPn_t));
                                sum2 = sum(chi(c).v_in_t(2:dat.BPn_t));
                                if sum2 > sum1,continue,end
                                sum1 = sum(chi(c+1).v_ex_t(2:dat.BPn_t));
                                sum2 = sum(chi(c+1).v_in_t(2:dat.BPn_t));
                                if sum2 > sum1,continue,end
                                
                                % Consertar o alinhamento do bordo de fuga (descomentar se o delta_z
                                % for usado como variável)
                                %chi(c).v_in(dat.BPn) = -chi(c).v_ex(dat.BPn)*chi(c).v_ex(dat.BPn-1)/chi(c).v_in(dat.BPn-1);
                                %chi(c+1).v_in(dat.BPn) = -chi(c+1).v_ex(dat.BPn)*chi(c+1).v_ex(dat.BPn-1)/chi(c+1).v_in(dat.BPn-1);
                            
                                chi(c).symm_t = 0;
                                chi(c+1).symm_t = 0;
                        
                            elseif chi(c).symm_t == 1 && chi(c+1).symm_t == 1 % Se ambos forem simétricos
                                n = randi([1,3]);
                                switch n                        
                                    case 1 % Trocar o raio do bordo de ataque
                                        
                                        temp1_ex = [af1_ex(1),af2_ex(2:end)];
                                        temp1_in = temp1_ex;
                                        temp2_ex = [af2_ex(1),af1_ex(2:end)];
                                        temp2_in = temp2_ex;
                                        
                                    case 2 % Trocar os pesos intermediários
                                        if dat.BPn_t == 2
                                            op = 1
                                        else
                                            op = randi([1,2]);
                                        end
                                        if op == 1 % Trocar tudo    
                                            temp1_ex = [af2_ex(1),af1_ex(2:dat.BPn_t),af2_ex(dat.BPn_t+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af1_ex(1),af2_ex(2:dat.BPn_t),af1_ex(dat.BPn_t+1:end)];
                                            temp2_in = temp2_ex;
                                            
                                        else % Trocar cortes
                                            num1 = randi(2:dat.BPn_t);
                                            temp1_1 = [af1_ex(2:num1),af2_ex(num1+1:dat.BPn_t)];
                                            temp2_1 = [af2_ex(2:num1),af1_ex(num1+1:dat.BPn_t)];
                                            temp1_ex = [af1_ex(1),temp2_1,af1_ex(dat.BPn_t+1:end)];
                                            temp1_in = temp1_ex;
                                            temp2_ex = [af2_ex(1),temp1_1,af2_ex(dat.BPn_t+1:end)];
                                            temp2_in = temp2_ex;

                                        end
                                        
                                    case 3 % Trocar os ângulos do bordo de fuga
                                        temp1_ex = [af2_ex(1:dat.BPn_t),af1_ex(dat.BPn_t+1),af2_ex(dat.BPn_t+2)];
                                        temp1_in = [af2_in(1:dat.BPn_t),af1_in(dat.BPn_t+1),af2_in(dat.BPn_t+2)];
                                        temp2_ex = [af1_ex(1:dat.BPn_t),af2_ex(dat.BPn_t+1),af1_ex(dat.BPn_t+2)];
                                        temp2_in = [af1_in(1:dat.BPn_t),af2_in(dat.BPn_t+1),af1_in(dat.BPn_t+2)];
                                        
                                end
                                
                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_t = temp1_ex;
                                    chi(c).v_in_t = temp1_in;
                                    chi(c+1).v_ex_t = temp2_ex;
                                    chi(c+1).v_in_t = temp2_in;
                                else
                                    chi(c).v_ex_t = temp2_ex;
                                    chi(c).v_in_t = temp2_in;
                                    chi(c+1).v_ex_t = temp1_ex;
                                    chi(c+1).v_in_t = temp1_in;
                                end
                                
                                chi(c).symm_t = 1;
                                chi(c+1).symm_t = 1;
                            
                            else % Se um for simétrico e o outro for assimétrico 
                                % Transformar o simétrico em um assimétrico e vice-versa
                                if chi(c).symm_t == 0 
                                    temp1_ex = af1_ex;
                                    temp1_in = af1_ex;
                                    temp2_ex = af2_ex;
                                    temp2_in = af1_in;
                                    
                                else
                                    temp1_ex = af2_ex;
                                    temp1_in = af2_ex;
                                    temp2_ex = af1_ex;
                                    temp2_in = af2_in;
                                    
                                end

                                % Decidir pra qual asa vai cada um dos aerofólios novos
                                if randi([0,1]) == 1 
                                    chi(c).v_ex_t = temp1_ex;
                                    chi(c).v_in_t = temp1_in;
                                    chi(c+1).v_ex_t = temp2_ex;
                                    chi(c+1).v_in_t = temp2_in;
                                else
                                    chi(c).v_ex_t = temp2_ex;
                                    chi(c).v_in_t = temp2_in;
                                    chi(c+1).v_ex_t = temp1_ex;
                                    chi(c+1).v_in_t = temp1_in;
                                end
                                
                                % Checar a separação dos bordos de fuga. Se não cumprirem o 
                                % requisito de separação, alterar o ângulo do bordo de fuga
                                % do intradorso
                                if chi(c).v_in_t(dat.BPn_t+1) < (dat.B_ext2_t(1)-chi(c).v_ex_t(dat.BPn_t+1))
                %                    chi(c).v_in(dat.BPn+1) = dat.B_ext2(1)-chi(c).v_ex(dat.BPn+1);
                                    chi(c).v_ex_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
                                    chi(c).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
                                end
                                if chi(c+1).v_in_t(dat.BPn_t+1) < (dat.B_ext2_t(1)-chi(c+1).v_ex_t(dat.BPn_t+1))
                                    chi(c+1).v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)-chi(c+1).v_ex_t(dat.BPn_t+1);
                                end
                                % Decisão de alterar o intradorso em base de uma nota na página
                                % 57(87) do Raymer (2018)
                                
                                % Checar os pesos
                                sum1 = sum(chi(c+1).v_ex_t(2:dat.BPn_t));
                                sum2 = sum(chi(c+1).v_in_t(2:dat.BPn_t));
                                if sum2 > sum1,continue,end
                                
                                chi(c).symm_t = 1;
                                chi(c+1).symm_t = 0;
                                
                            end
                    end 
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(chi(c).v_ex_t,chi(c).v_in_t,dat,[dat.N1_t,dat.N2_t]),dat);
                    if check == 0,continue,end
                    check = quality(run_cst_TCC2_3D(chi(c+1).v_ex_t,chi(c+1).v_in_t,dat,[dat.N1_t,dat.N2_t]),dat); 
                
                end
            end
            
            % Torção geométrica na ponta
            if cross_op_v(9) == 1 
                chi(c).tw_t = pop(par(2)).tw_t;
                chi(c+1).tw_t = pop(par(1)).tw_t;
            end
            
            % Enflechamento total (asas trapezoidais simples)
            if cross_op_v(10) == 1 && pop(par(1)).type == 0 && pop(par(2)).type == 0
                chi(c).sweep = pop(par(2)).sweep;
                chi(c+1).sweep = pop(par(1)).sweep;
            end 
        
            % Enflechamento da primeira seção (asas bitrapezoidais)
            if cross_op_v(11) == 1 && pop(par(1)).type == 1 && pop(par(2)).type == 1
                chi(c).sweep1 = pop(par(2)).sweep1;
                chi(c+1).sweep1 = pop(par(1)).sweep1;
            end
            
            % Enflechamento da segunda seção (asas bitrapezoidais)
            if cross_op_v(12) == 1 && pop(par(1)).type == 1 && pop(par(2)).type == 1
                chi(c).sweep2 = pop(par(2)).sweep2;
                chi(c+1).sweep2 = pop(par(1)).sweep2;
            end
        
%        else % Transformar uma asa trapezoidal simples em trapezoidal dupla e vice-versa
%             % (adiciona-se/retira-se o aerofólio do meio)
%            if pop(par(1)).type == 0 % Se a primeira do par for trapezoidal simples
%                chi(c) = pop(par(1));
%                chi(c).type = 1;
%                chi(c).b1 = pop(par(2)).b1;
%                chi(c).c_m = pop(par(2)).c_m;
%                chi(c).v_ex_m = pop(par(2)).v_ex_m;
%                chi(c).v_in_m = pop(par(2)).v_in_m;
%                chi(c).tw_m = pop(par(2)).tw_m;
%                chi(c).sweep1 = pop(par(2)).sweep1;	
%                chi(c).sweep2 = pop(par(2)).sweep2;
%                
%                chi(c+1) = pop(par(2));
%                chi(c+1).type = 0;
%                chi(c+1).b1 = [];
%                chi(c+1).c_m = [];
%                chi(c+1).v_ex_m = zeros(1,dat.BPn_m+2);
%                chi(c+1).v_in_m = zeros(1,dat.BPn_m+2);
%                chi(c+1).tw_m = [];
%                chi(c+1).sweep = pop(par(1)).sweep;
%                
%                
%            else % Se a primeira do par for trapezoidal dupla
%                chi(c) = pop(par(2));
%                chi(c).type = 1;
%                chi(c).b1 = pop(par(1)).b1;
%                chi(c).c_m = pop(par(1)).c_m;
%                chi(c).v_ex_m = pop(par(1)).v_ex_m;
%                chi(c).v_in_m = pop(par(1)).v_in_m;
%%                chi(c).af_m = pop(par(1)).af_m;
%                chi(c).tw_m = pop(par(1)).tw_m;
%                chi(c).sweep1 = pop(par(1)).sweep1;	
%                chi(c).sweep2 = pop(par(1)).sweep2;
%                
%                chi(c+1) = pop(par(1));
%                chi(c+1).type = 0;
%                chi(c+1).b1 = [];
%                chi(c+1).c_m = [];
%                chi(c+1).v_ex_m = zeros(1,dat.BPn_m+2);
%                chi(c+1).v_in_m = zeros(1,dat.BPn_m+2);
%                chi(c+1).tw_m = [];
%                chi(c+1).sweep = pop(par(2)).sweep;
%                
%            end
%            
%            % Checar o requisito das envergaduras
%            chi(c).b1 = dat.or_b1 + rand*dat.b1_ext(randi(2));
%            if chi(c).b1 < dat.b1_ext(3)
%                chi(c).b1 = dat.b1_ext(3);
%            end
%            if chi(c).b1 > chi(c).b - 2*dat.b1_ext(4)
%                chi(c).b1 = chi(c).b - 2*dat.b1_ext(4);
%            end
%            
%            % Checar o requisito dos comprimentos de corda
%            if ~dat.or_c_m == 'L'
%                if chi(c).c_t > chi(c).c_m || chi(c).c_m > chi(c).c_r
%                    chi(c).c_m = dat.or_c_m + rand*dat.c_m_ext(randi(2));
%                    if chi(c).c_m < dat.c_m_ext(3)
%                        chi(c).c_m = dat.c_m_ext(3);
%                    end
%                    if chi(c).c_m > chi(c).c_r
%                        chi(c).c_m = chi(c).c_r;
%                    end
%                    if chi(c).c_t > chi(c).c_m
%                        chi(c).c_t = chi(c).c_m;
%                    end
%                end
%            end
            
%            if chi(c).c_t > chi(c).c_m || chi(c).c_m > chi(c).c_r
%                % Se não cumprir o requisito, atribuir um novo valor
%                dat.c_m_ext = [0,0];
%                if chi(c).c_t >= dat.c_m_ext_in(1)
%                    dat.c_m_ext(1) = chi(c).c_t;
%                else
%                    dat.c_m_ext(1) = dat.c_m_ext_in(1);
%                end
%                if chi(c).c_r <= dat.c_m_ext_in(2)
%                    dat.c_m_ext(2) = chi(c).c_r;
%                else
%                    dat.c_m_ext(2) = dat.c_m_ext_in(2);
%                end
%                dat.c_m_ext = dat.c_m_ext(1):dat.c_m_step:dat.c_m_ext(2);
%                chi(c).c_m = dat.c_m_ext(randi(length(dat.c_m_ext)));
%            end
            
        end
        
        % Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
        if chi(c).type == 1
            if chi(c).b1 > chi(c).b-2*dat.b1_ext(4) || chi(c).c_t > chi(c).c_m && dat.or_c_m ~= 'L' || chi(c).c_m > chi(c).c_r && dat.or_c_m ~= 'L' 
                error('Problema na geometria da planta (c)')
            end 
        end
        if chi(c+1).type == 1
            if chi(c+1).b1 > chi(c+1).b-2*dat.b1_ext(4) || chi(c+1).c_t > chi(c+1).c_m && dat.or_c_m ~= 'L' || chi(c+1).c_m > chi(c+1).c_r && dat.or_c_m ~= 'L' 
                error('Problema na geometria da planta (c+1)')
            end 
        end
        
        % Debugging: ver se aerofólios assimétricos são transformados em simétricos por algum motivo
        if chi(c).symm_r ~= 1 && chi(c).v_ex_r == chi(c).v_in_r
            error('Aerofólio assimétrico transformado em simétrico (raiz)')
        end
        if chi(c).type == 1 && chi(c).symm_m ~= 1 && chi(c).v_ex_m == chi(c).v_in_m
            error('Aerofólio assimétrico transformado em simétrico (meio)')
        end
        if chi(c).symm_t ~= 1 && chi(c).v_ex_t == chi(c).v_in_t
            error('Aerofólio assimétrico transformado em simétrico (ponta)')
        end
        if chi(c+1).symm_r ~= 1 && chi(c+1).v_ex_r == chi(c+1).v_in_r
            error('Aerofólio assimétrico transformado em simétrico (raiz)')
        end
        if chi(c+1).type == 1 && chi(c+1).symm_m ~= 1 && chi(c+1).v_ex_m == chi(c+1).v_in_m
            error('Aerofólio assimétrico transformado em simétrico (meio)')
        end
        if chi(c+1).symm_t ~= 1 && chi(c+1).v_ex_t == chi(c+1).v_in_t
            error('Aerofólio assimétrico transformado em simétrico (ponta)')
        end
        
        % Debugging: ver se a opção 'L' da corda do meio está sendo cumprida
        if dat.or_c_m == 'L' && chi(c).c_m ~= 'L'
            error("Erro da opção 'L' da corda do meio (c)")
        end
        if dat.or_c_m == 'L' && chi(c+1).c_m ~= 'L'
            error("Erro da opção 'L' da corda do meio (c+1)")
        end
        
        % Zerar a pontuação dos filhos 
        chi(c).score = 0;
        chi(c+1).score = 0;
        
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
            
            
            % A mutação funciona de modo a atribuir novos valores aos respectivos campos
            
            % Genes a serem alterados 
            % Envergadura b
            % Envergadura b1 (asas bitrapezoidais apenas)
            % Corda da raiz
            % Corda do meio (asas bitrapezoidais apenas)
            % Corda da ponta
            % - aerofólio da raiz
            % - aerofólio do meio (asas bitrapezoidais apenas)
            % - aerofólio da ponta
            % Torção geométrica na ponta
            
            mu_op_v = randi([0,1],1,12);
            
            % Envergadura b
            if mu_op_v(1) == 1
                chi(s).b = dat.or_b + rand*dat.b_ext(randi(2));
                if chi(s).b < dat.b_ext(3)
                    chi(s).b = dat.b_ext(3);
                end
            end
            
            % Envergadura da primeira seção b1 (asas bitrapezoidais apenas)
            if mu_op_v(2) == 1 && chi(s).type == 1
                chi(s).b1 = dat.or_b1 + rand*dat.b1_ext(randi(2));
                if chi(s).b1 < dat.b1_ext(3)
                    chi(s).b1 = dat.b1_ext(3);
                end
            end
            
            % Corda da raiz c_r
            if mu_op_v(3) == 1
                chi(s).c_r = dat.or_c_r + rand*dat.c_r_ext(randi(2));
                if chi(s).c_r < dat.c_r_ext(3)
                    chi(s).c_r = dat.c_r_ext(3);
                end
				
				% Se a corda da ponta for maior do que a da raiz, atribuir o 
                % valor da ponta à raiz
                if chi(s).c_t > chi(s).c_r
                    chi(s).c_r = chi(s).c_t;
                end
				
            end
            
            % Corda da ponta c_t (deve cumprir o requisito c_t <= c_r)
            if mu_op_v(4) == 1
                chi(s).c_t = dat.or_c_t + rand*dat.c_t_ext(randi(2));
                if chi(s).c_t < dat.c_t_ext(3)
                    chi(s).c_t = dat.c_t_ext(3);
                end
				
				% Se a corda da ponta for maior do que a da raiz, atribuir o 
                % valor da raiz à ponta
                if chi(s).c_t > chi(s).c_r
                    chi(s).c_t = chi(s).c_r;
                end
				
            end
            
            % Corda do meio c_m (asas bitrapezoidais apenas)
            if mu_op_v(5) == 1 && ~dat.or_c_m == 'L' 
                chi(s).c_m = dat.or_c_m + rand*dat.c_m_ext(randi(2));
                if chi(s).c_m < dat.c_m_ext(3)
                    chi(s).c_m = dat.c_m_ext(3);
                end
            end
            
            % Correções referentes à planta da asa (asas bitrapezoidais apenas)
            if chi(s).type == 1
                
                % Se a envergadura da primeira seção for maior do que permitido
                % pelo requisito de separação, atribuir o máximo valor que 
                % cumpre o requisito
                if chi(s).b1 > chi(s).b-dat.b1_ext(4)*2
                    chi(s).b1 = chi(s).b-dat.b1_ext(4)*2;
                end
                
                % Se a corda do meio for maior do que a corda da raiz, atribuir
                % o valor da raiz ao meio
                if chi(s).c_m > chi(s).c_r && dat.or_c_m ~= 'L'
                    chi(s).c_m = chi(s).c_r;
                end
                
                % Se a corda do meio for menor que a corda da ponta, atribuir
                % o valor da ponta ao meio
                if chi(s).c_m < chi(s).c_t && dat.or_c_m ~= 'L'
                    chi(s).c_m = chi(s).c_t;
                end
                
            end
            
            % Aerofólio da raiz
            if mu_op_v(6) == 1
                check = 0;
                while check == 0
                    
                    temp = chi(s);
                    if temp.symm_r == 1 % Caso o perfil seja simétrico
                        n = [1,4,5](randi(3));
                    else % Caso o perfil seja assimétrico
                        n = randi(5);
                    end
                    switch n
                        case 1 % Alterar o raio do bordo de ataque
                            if temp.symm_r == 1
                                P = 1;
                            else
                                P = randi([1,4]);
                            end
                            if P == 1 % Mudar ambos para o mesmo valor
                                temp.v_ex_r(1) = dat.or_v_ex_r(1) + rand*dat.le_R_ext1_r(randi(2));
                                temp.v_in_r(1) = temp.v_ex_r(1);
                            elseif P == 2 % Mudar ambos independentemente 
                                temp.v_ex_r(1) = dat.or_v_ex_r(1) + rand*dat.le_R_ext1_r(randi(2));
                                temp.v_in_r(1) = dat.or_v_in_r(1) + rand*dat.le_R_ext2_r(randi(2));
                            elseif P == 3 % Mudar do extradorso
                                temp.v_ex_r(1) = dat.or_v_ex_r(1) + rand*dat.le_R_ext1_r(randi(2));
                            else % Mudar do intradorso
                                temp.v_in_r(1) = dat.or_v_in_r(1) + rand*dat.le_R_ext2_r(randi(2));
                            end
							
							% Checar o raio do bordo de ataque (simétricos)
							if temp.symm_r == 1 && temp.v_ex_r(1) < dat.le_R_ext1_r(3)
								temp.v_ex_r(1) = dat.le_R_ext1_r(3);
								temp.v_in_r(1) = temp.v_ex_r(1);
							end
							
							% Checar o raio do bordo de ataque (assimétricos)
							if temp.symm_r == 0 && temp.v_ex_r(1) < dat.le_R_ext1_r(3)
								temp.v_ex_r(1) = dat.le_R_ext1_r(3);
							end
							if temp.symm_r == 0 && temp.v_in_r(1) < dat.le_R_ext2_r(3)
								temp.v_in_r(1) = dat.le_R_ext2_r(3);
							end
                            
                        case 2 % Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_r
                                temp.v_ex_r(a) = dat.or_v_ex_r(a) + rand*dat.A_ext1_r(randi(2));
                            end
                            
                        case 3 % Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_r
                                temp.v_in_r(a) = dat.or_v_in_r(a) + rand*dat.A_ext2_r(randi(2));
                            end
                            
                        case 4 % Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                            if temp.symm_r == 1 % Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                                num = rand(1,dat.BPn_r-1)*0.1;
                                for a = 2:dat.BPn_r
                                    temp.v_ex_r(a) = dat.or_v_ex_r(a) + rand*dat.A_ext1_r(randi(2));
                                    temp.v_in_r(a) = temp.v_ex_r(a);
                                end
                            else % Perfis assimétricos ficam com pesos intermediários distintos
                                for a = 2:dat.BPn_r
                                    temp.v_ex_r(a) = dat.or_v_ex_r(a) + rand*dat.A_ext1_r(randi(2));
                                end
                                
                                for a = 2:dat.BPn_r
                                    temp.v_in_r(a) = dat.or_v_in_r(a) + rand*dat.A_ext2_r(randi(2));
                                end
                            end
                            
                        case 5 % Alterar o ângulo do bordo de fuga
%                            if sum(dat.B_ext1_r) ~= 0 && dat.B_ext2_r(2) ~= 0
                            if temp.symm_r == 1 % Se for simétrico
								temp.v_ex_r(dat.BPn_r+1) = dat.or_v_ex_r(dat.BPn_r+1) + rand*dat.B_ext1_r(randi(2));
								temp.v_in_r(dat.BPn_r+1) = temp.v_ex_r(dat.BPn_r+1);
								
								% Checar separação do bordo de fuga
								if temp.v_ex_r(dat.BPn_r+1) + temp.v_in_r(dat.BPn_r+1) < dat.B_ext2_r(1)
									temp.v_ex_r(dat.BPn_r+1) = dat.B_ext2_r(1)/2;
									temp.v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1)/2;
								end
							else % Se for assimétrico
                                % Ignorar a mutação de beta pra perfis assimétricos caso
                                % a intenção seja não alterar o ângulo de cima. Devido à
                                % forma que o ângulo do intradorso foi implementado, ele
                                % necessariamente depende do ângulo do extradorso. O que 
                                % se deve fazer depois é permitir a manipulação do ângulo
                                % do intradorso igual ao que é feito no extradorso, mas 
                                % em seguida fazendo checagens e alterando o valor se 
                                % preciso (como os requisitos das plantas das asas).
                                % No caso do não cumprimento do requisito, deve-se alterar
                                % ou apenas um dos ângulos ou ambos igualmente.
                                if sum(dat.B_ext1_r) ~= 0 
                                    
                                    temp.v_ex_r(dat.BPn_r+1) = dat.or_v_ex_r(dat.BPn_r+1) + rand*dat.B_ext1_r(randi(2));
                                    temp.v_in_r(dat.BPn_r+1) = (dat.B_ext2_r(1) - temp.v_ex_r(dat.BPn_r+1)) + rand*dat.B_ext2_r(2);
                                    
                                    % Checar separação do bordo de fuga (beta2>=L-beta1)
                                    if temp.v_in_r(dat.BPn_r+1) < dat.B_ext2_r(1) - temp.v_ex_r(dat.BPn_r+1)
                                        temp.v_in_r(dat.BPn_r+1) = dat.B_ext2_r(1) - temp.v_ex_r(dat.BPn_r+1);
                                    end
                                end
							end
							
                    end
                     
                    % Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    % igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_r(2:dat.BPn_r));
                    sum2 = sum(temp.v_in_r(2:dat.BPn_r));
                    if sum2 > sum1
                        continue
                    end
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_r,temp.v_in_r,dat,[dat.N1_r,dat.N2_r]),dat); 
                end
                chi(s) = temp;
                
            end
            
            % Aerofólio do meio (asas bitrapezoidais apenas)
            if mu_op_v(7) == 1 && chi(s).type == 1 && ~ismember('L',dat.le_R_ext1_m)
                check = 0;
                while check == 0
                    
                    temp = chi(s);
                    if temp.symm_m == 1 % Caso o perfil seja simétrico
                        n = [1,4,5](randi(3));
                    else % Caso o perfil seja assimétrico
                        n = randi(5);
                    end
                    switch n
                        case 1 % Alterar o raio do bordo de ataque
                            if temp.symm_m == 1
                                P = 1;
                            else
                                P = randi([1,4]);
                            end
                            if P == 1 % Mudar ambos para o mesmo valor
                                temp.v_ex_m(1) = dat.or_v_ex_m(1) + rand*dat.le_R_ext1_m(randi(2));
                                temp.v_in_m(1) = temp.v_ex_m(1);
                            elseif P == 2 % Mudar ambos independentemente 
                                temp.v_ex_m(1) = dat.or_v_ex_m(1) + rand*dat.le_R_ext1_m(randi(2));
                                temp.v_in_m(1) = dat.or_v_in_m(1) + rand*dat.le_R_ext2_m(randi(2));
                            elseif P == 3 % Mudar do extradorso
                                temp.v_ex_m(1) = dat.or_v_ex_m(1) + rand*dat.le_R_ext1_m(randi(2));
                            else % Mudar do intradorso
                                temp.v_in_m(1) = dat.or_v_in_m(1) + rand*dat.le_R_ext2_m(randi(2));
                            end
							
							% Checar o raio do bordo de ataque (simétricos)
							if temp.symm_m == 1 && temp.v_ex_m(1) < dat.le_R_ext1_m(3)
								temp.v_ex_m(1) = dat.le_R_ext1_m(3);
								temp.v_in_m(1) = temp.v_ex_m(1);
							end
							
							% Checar o raio do bordo de ataque (assimétricos)
							if temp.symm_m == 0 && temp.v_ex_m(1) < dat.le_R_ext1_m(3)
								temp.v_ex_m(1) = dat.le_R_ext1_m(3);
							end
							if temp.symm_m == 0 && temp.v_in_m(1) < dat.le_R_ext2_m(3)
								temp.v_in_m(1) = dat.le_R_ext2_m(3);
							end
                            
                        case 2 % Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_m
                                temp.v_ex_m(a) = dat.or_v_ex_m(a) + rand*dat.A_ext1_m(randi(2));
                            end
                            
                        case 3 % Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_m
                                temp.v_in_m(a) = dat.or_v_in_m(a) + rand*dat.A_ext2_m(randi(2));
                            end
                            
                        case 4 % Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                            if temp.symm_m == 1 % Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                                num = rand(1,dat.BPn_m-1)*0.1;
                                for a = 2:dat.BPn_m
                                    temp.v_ex_m(a) = dat.or_v_ex_m(a) + rand*dat.A_ext1_m(randi(2));
                                    temp.v_in_m(a) = temp.v_ex_m(a);
                                end
                            else % Perfis assimétricos ficam com pesos intermediários distintos
                                for a = 2:dat.BPn_m
                                    temp.v_ex_m(a) = dat.or_v_ex_m(a) + rand*dat.A_ext1_m(randi(2));
                                end
                                
                                for a = 2:dat.BPn_m
                                    temp.v_in_m(a) = dat.or_v_in_m(a) + rand*dat.A_ext2_m(randi(2));
                                end
                            end
                            
                        case 5 % Alterar o ângulo do bordo de fuga
                            if temp.symm_m == 1 % Se for simétrico
								temp.v_ex_m(dat.BPn_m+1) = dat.or_v_ex_m(dat.BPn_m+1) + rand*dat.B_ext1_m(randi(2));
								temp.v_in_m(dat.BPn_m+1) = temp.v_ex_m(dat.BPn_m+1);
								
								% Checar separação do bordo de fuga
								if temp.v_ex_m(dat.BPn_m+1) + temp.v_in_m(dat.BPn_m+1) < dat.B_ext2_m(1)
									temp.v_ex_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
									temp.v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1)/2;
								end
							else % Se for assimétrico
								temp.v_ex_m(dat.BPn_m+1) = dat.or_v_ex_m(dat.BPn_m+1) + rand*dat.B_ext1_m(randi(2));
								temp.v_in_m(dat.BPn_m+1) = (dat.B_ext2_m(1) - temp.v_ex_m(dat.BPn_m+1)) + rand*dat.B_ext2_m(2);
								
								% Checar separação do bordo de fuga (beta2>=L-beta1)
								if temp.v_in_m(dat.BPn_m+1) < dat.B_ext2_m(1) - temp.v_ex_m(dat.BPn_m+1)
									temp.v_in_m(dat.BPn_m+1) = dat.B_ext2_m(1) - temp.v_ex_m(dat.BPn_m+1);
								end
							end
							
                    end
                     
                    % Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    % igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_m(2:dat.BPn_m));
                    sum2 = sum(temp.v_in_m(2:dat.BPn_m));
                    if sum2 > sum1
                        continue
                    end
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_m,temp.v_in_m,dat,[dat.N1_m,dat.N2_m]),dat); 
                end
                chi(s) = temp;
            end
            
            % Aerofólio da ponta
            if mu_op_v(8) == 1
                check = 0;
                while check == 0
                    
                    temp = chi(s);
                    if temp.symm_t == 1 % Caso o perfil seja simétrico
                        n = [1,4,5](randi(3));
                    else % Caso o perfil seja assimétrico
                        n = randi(5);
                    end
                    switch n
                        case 1 % Alterar o raio do bordo de ataque
                            if temp.symm_t == 1
                                P = 1;
                            else
                                P = randi([1,4]);
                            end
                            if P == 1 % Mudar ambos para o mesmo valor
                                temp.v_ex_t(1) = dat.or_v_ex_t(1) + rand*dat.le_R_ext1_t(randi(2));
                                temp.v_in_t(1) = temp.v_ex_t(1);
                            elseif P == 2 % Mudar ambos independentemente 
                                temp.v_ex_t(1) = dat.or_v_ex_t(1) + rand*dat.le_R_ext1_t(randi(2));
                                temp.v_in_t(1) = dat.or_v_in_t(1) + rand*dat.le_R_ext2_t(randi(2));
                            elseif P == 3 % Mudar do extradorso
                                temp.v_ex_t(1) = dat.or_v_ex_t(1) + rand*dat.le_R_ext1_t(randi(2));
                            else % Mudar do intradorso
                                temp.v_in_t(1) = dat.or_v_in_t(1) + rand*dat.le_R_ext2_t(randi(2));
                            end
							
							% Checar o raio do bordo de ataque (simétricos)
							if temp.symm_t == 1 && temp.v_ex_t(1) < dat.le_R_ext1_t(3)
								temp.v_ex_t(1) = dat.le_R_ext1_t(3);
								temp.v_in_t(1) = temp.v_ex_t(1);
							end
							
							% Checar o raio do bordo de ataque (assimétricos)
							if temp.symm_t == 0 && temp.v_ex_t(1) < dat.le_R_ext1_t(3)
								temp.v_ex_t(1) = dat.le_R_ext1_t(3);
							end
							if temp.symm_t == 0 && temp.v_in_t(1) < dat.le_R_ext2_t(3)
								temp.v_in_t(1) = dat.le_R_ext2_t(3);
							end
                            
                        case 2 % Alterar os pesos intermediários (extradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_t
                                temp.v_ex_t(a) = dat.or_v_ex_t(a) + rand*dat.A_ext1_t(randi(2));
                            end
                            
                        case 3 % Alterar os pesos intermediários (intradorso) dentro de uma extensão próxima aos valores originais
                            for a = 2:dat.BPn_t
                                temp.v_in_t(a) = dat.or_v_in_t(a) + rand*dat.A_ext2_t(randi(2));
                            end
                            
                        case 4 % Alterar os pesos intermediários (extradorso e intradorso) dentro de uma extensão próxima aos valores originais
                            if temp.symm_t == 1 % Perfis simétricos ficam com os pesos intermediários com os mesmos valores
                                num = rand(1,dat.BPn_t-1)*0.1;
                                for a = 2:dat.BPn_t
                                    temp.v_ex_t(a) = dat.or_v_ex_t(a) + rand*dat.A_ext1_t(randi(2));
                                    temp.v_in_t(a) = temp.v_ex_t(a);
                                end
                            else % Perfis assimétricos ficam com pesos intermediários distintos
                                for a = 2:dat.BPn_t
                                    temp.v_ex_t(a) = dat.or_v_ex_t(a) + rand*dat.A_ext1_t(randi(2));
                                end
                                
                                for a = 2:dat.BPn_t
                                    temp.v_in_t(a) = dat.or_v_in_t(a) + rand*dat.A_ext2_t(randi(2));
                                end
                            end
                            
                        case 5 % Alterar o ângulo do bordo de fuga
							if temp.symm_t == 1 % Se for simétrico
								temp.v_ex_t(dat.BPn_t+1) = dat.or_v_ex_t(dat.BPn_t+1) + rand*dat.B_ext1_t(randi(2));
								temp.v_in_t(dat.BPn_t+1) = temp.v_ex_t(dat.BPn_t+1);
								
								% Checar separação do bordo de fuga
								if temp.v_ex_t(dat.BPn_t+1) + temp.v_in_t(dat.BPn_t+1) < dat.B_ext2_t(1)
									temp.v_ex_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
									temp.v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1)/2;
								end
							else % Se for assimétrico
								temp.v_ex_t(dat.BPn_t+1) = dat.or_v_ex_t(dat.BPn_t+1) + rand*dat.B_ext1_t(randi(2));
								temp.v_in_t(dat.BPn_t+1) = (dat.B_ext2_t(1) - temp.v_ex_t(dat.BPn_t+1)) + rand*dat.B_ext2_t(2);
								
								% Checar separação do bordo de fuga (beta2>=L-beta1)
								if temp.v_in_t(dat.BPn_t+1) < dat.B_ext2_t(1) - temp.v_ex_t(dat.BPn_t+1)
									temp.v_in_t(dat.BPn_t+1) = dat.B_ext2_t(1) - temp.v_ex_t(dat.BPn_t+1);
								end
							end

                    end
                     
                    % Checar os pesos (soma de pesos do intradorso deve ser menor ou
                    % igual à soma de pesos do extradorso)
                    sum1 = sum(temp.v_ex_t(2:dat.BPn_t));
                    sum2 = sum(temp.v_in_t(2:dat.BPn_t));
                    if sum2 > sum1
                        continue
                    end
                    
                    % Checagem de qualidade
                    check = quality(run_cst_TCC2_3D(temp.v_ex_t,temp.v_in_t,dat,[dat.N1_t,dat.N2_t]),dat); 
                end
                chi(s) = temp;
            end
            
            % Torção geométrica na ponta
            if mu_op_v(9) == 1
                chi(s).tw_t = dat.or_tw_t + rand*dat.tw_t_ext(randi(2));
            end
            
            % Enflechamento total (asas trapezoidais simples)
            if mu_op_v(10) == 1 && chi(s).type == 0 && ~ismember('Z',dat.sweep_ext) && dat.or_sweep ~= 'Z'
                chi(s).sweep = dat.or_sweep + rand*dat.sweep_ext(randi(2));
            end
            
            % Enflechamento da primeira seção (asas trapezoidais duplas)
            if mu_op_v(11) == 1 && chi(s).type == 1 && ~ismember('Z',dat.sweep1_ext) && dat.or_sweep1 ~= 'Z'
                chi(s).sweep1 = dat.or_sweep1 + rand*dat.sweep1_ext(randi(2));
            end
            
            % Enflechamento da segunda seção (asas trapezoidais duplas)
            if mu_op_v(12) == 1 && chi(s).type == 1 && ~ismember('Z',dat.sweep2_ext) && dat.or_sweep2 ~= 'Z'
                chi(s).sweep2 = dat.or_sweep2 + rand*dat.sweep2_ext(randi(2));
            end
            
            % Debugging: ver como estão sendo cumpridas os requisitos de geometria da planta
			if chi(s).type == 1
				if chi(s).b1 > chi(s).b-2*dat.b1_ext(4) || chi(s).c_t > chi(s).c_m && dat.or_c_m~='L' || chi(s).c_m > chi(s).c_r && dat.or_c_m~='L'
					error('Problema na geometria da planta')
				end
			end
        
            % Debugging: ver se aerofólios assimétricos são transformados em simétricos por algum motivo
            if chi(s).symm_r ~= 1 && chi(s).v_ex_r == chi(s).v_in_r
                error('Aerofólio assimétrico transformado em simétrico (raiz)')
            end
            if chi(s).type == 1 && chi(s).symm_m ~= 1 && chi(s).v_ex_m == chi(s).v_in_m
                error('Aerofólio assimétrico transformado em simétrico (meio)')
            end
            if chi(s).symm_t ~= 1 && chi(s).v_ex_t == chi(s).v_in_t
                error('Aerofólio assimétrico transformado em simétrico (ponta)')
            end
            
            if dat.or_c_m == 'L' && chi(s).c_m ~= 'L'
                error("Erro da opção 'L' da corda do meio")
            end
            
            k = k+1;
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

clc
t = toc;
min = fix(t/60); s = rem(t,60);
fprintf('Tempo: %d min e %.2f s\n\n', min, s)

% Comparar aerofólios
figure(2),clf,hold on
%plot_airfoil_cst_TCC2(run_cst_TCC2(pop(pos).v_ex,pop(pos).v_in,dat),1,loop),hold on
plot_airfoil_cst_TCC2(run_cst_TCC2_3D(pop(pos).v_ex_r,pop(pos).v_in_r,dat,[dat.N1_r,dat.N2_r]),4)
plot_airfoil_cst_TCC2(run_cst_TCC2_3D(dat.or_v_ex_r,dat.or_v_in_r,dat,[dat.N1_r,dat.N2_r]),3)
legend('Otimizado','Original'),title('Raiz'),axis equal,grid on
if dat.type == 1
    figure(3),clf,hold on
    plot_airfoil_cst_TCC2(run_cst_TCC2_3D(pop(pos).v_ex_m,pop(pos).v_in_m,dat,[dat.N1_m,dat.N2_m]),4)
    plot_airfoil_cst_TCC2(run_cst_TCC2_3D(dat.or_v_ex_m,dat.or_v_in_m,dat,[dat.N1_m,dat.N2_m]),3)
    legend('Otimizado','Original'),title('Meio'),axis equal,grid on  
end
figure(4),clf,hold on
plot_airfoil_cst_TCC2(run_cst_TCC2_3D(pop(pos).v_ex_t,pop(pos).v_in_t,dat,[dat.N1_t,dat.N2_t]),4)
plot_airfoil_cst_TCC2(run_cst_TCC2_3D(dat.or_v_ex_t,dat.or_v_in_t,dat,[dat.N1_r,dat.N2_t]),3)
legend('Otimizado','Original'),title('Ponta'),axis equal,grid on

% Fazer gráficos dos coeficientes dos melhores indivíduos de cada geração
for i = 1:dat.cases
    aero_m = zeros(dat.iter,4);
    for j = 1:dat.iter
        aero_m(j,:) = archive(j).aero(i,:);
    end
    figure(i+4),clf,hold on,grid on
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
    
    % Trocar separador decimal
    xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
    new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
    set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

end

% Pegar o struct de arquivo e imprimir todos
for i = 1:length(archive)
    disp(['<< Iteração ' num2str(i) ' >>'])
    disp(' - Dados da planta - ')
    if archive(i).type == 0
        disp('Tipo: trapezoidal simples (type = 0)')
    else
        disp('Tipo: trapezoidal dupla (type = 1)')
    end
    fprintf('Envergadura b = %.6f\n',archive(i).b)
    if archive(i).type == 1,fprintf('Envergadura da primeira seção b1 = %.6f\n',archive(i).b1),end
    fprintf('Corda da raiz c_r = %.6f\n',archive(i).c_r)
    if archive(i).type == 1,fprintf('Corda do meio c_m = %.6f\n',archive(i).c_m),end
    fprintf('Corda da ponta c_t = %.6f\n\n',archive(i).c_t)
    if archive(i).type == 0
        fprintf('Enflechamento completo sweep = %.6f\n\n',archive(i).sweep)
    else
        fprintf('Enflechamento da primeira seção sweep1 = %.6f°\n',archive(i).sweep1)
        fprintf('Enflechamento da segunda seção sweep2 = %.6f°\n\n',archive(i).sweep2)
    end
       
    
    disp('- Dados dos aerofólios -')
    
    disp(['Raiz (Polinômio de grau ' num2str(dat.BPn_r) '):'])
    fprintf('v_ex = [%.4f, ', archive(i).v_ex_r(1))
    for j = 2:(length(pop(1).v_ex_r)-2)
        fprintf('%.4f, ',archive(i).v_ex_r(j))
    end
    fprintf('%.4f, ', archive(i).v_ex_r(end-1))
    fprintf('%.4f];\n', archive(i).v_ex_r(end))
    
    fprintf('v_in = [%.4f, ', archive(i).v_in_r(1))
    for j = 2:(length(pop(1).v_in_r)-2)
        fprintf('%.4f, ',archive(i).v_in_r(j))
    end
    fprintf('%.4f, ', archive(i).v_in_r(end-1))
    fprintf('%.4f];\n\n', archive(i).v_in_r(end))
    
    if archive(i).type == 1
        disp(['Meio (Polinômio de grau ' num2str(dat.BPn_m) '):'])
        fprintf('v_ex = [%.4f, ', archive(i).v_ex_m(1))
        for j = 2:(length(pop(1).v_ex_m)-2)
            fprintf('%.4f, ',archive(i).v_ex_m(j))
        end
        fprintf('%.4f, ', archive(i).v_ex_m(end-1))
        fprintf('%.4f];\n', archive(i).v_ex_m(end))
        
        fprintf('v_in = [%.4f, ', archive(i).v_in_m(1))
        for j = 2:(length(pop(1).v_in_m)-2)
            fprintf('%.4f, ',archive(i).v_in_m(j))
        end
        fprintf('%.4f, ', archive(i).v_in_m(end-1))
        fprintf('%.4f];\n\n', archive(i).v_in_m(end))
    end
    
    disp(['Ponta (Polinômio de grau ' num2str(dat.BPn_t) '):'])
    fprintf('v_ex = [%.4f, ', archive(i).v_ex_t(1))
    for j = 2:(length(pop(1).v_ex_t)-2)
        fprintf('%.4f, ',archive(i).v_ex_t(j))
    end
    fprintf('%.4f, ', archive(i).v_ex_t(end-1))
    fprintf('%.4f];\n', archive(i).v_ex_t(end))
    
    fprintf('v_in = [%.4f, ', archive(i).v_in_t(1))
    for j = 2:(length(pop(1).v_in_t)-2)
        fprintf('%.4f, ',archive(i).v_in_t(j))
    end
    fprintf('%.4f, ', archive(i).v_in_t(end-1))
    fprintf('%.4f];\n\n', archive(i).v_in_t(end))
            
    fprintf('Torção geométrica na ponta tw_t = %.6f\n\n',archive(i).tw_t)
    
    disp('- Dados aerodinâmicos -')
    for j = 1:dat.cases
%        fprintf('Condição de voo %d (CL,CD,L/D,CM): ',j),disp(archive(i).aero(j,:))
        fprintf('Condição de voo %d (CL,CD,L/D,CM)\nOriginal: ',j),disp(original.aero(j,:))
        fprintf('Otimizado: '),disp(archive(i).aero(j,:))
        
        if ismember('q',dat.coeff_op(:,1)) || ismember('#',dat.coeff_op(:,1))
            fprintf('L = %f N (Original)\n',original.aero(j,1)*1/2*dat.rho(j)*dat.v_ref(j)^2*original(i).S)
            fprintf('L = %f N (Otimizado)\n',archive(i).aero(j,1)*1/2*dat.rho(j)*dat.v_ref(j)^2*archive(i).S)
        end
        
        if ismember('q',dat.coeff_op(:,2)) || ismember('#',dat.coeff_op(:,2))
            fprintf('D = %f N (Original)\n',original.aero(j,2)*1/2*dat.rho(j)*dat.v_ref(j)^2*original(i).S)
            fprintf('D = %f N (Otimizado)\n',archive(i).aero(j,2)*1/2*dat.rho(j)*dat.v_ref(j)^2*archive(i).S)
        end
        
        if ismember('q',dat.coeff_op(:,4)) 
            fprintf('M = %f Nm (Original)\n',original.aero(j,4)*1/2*dat.rho(j)*dat.v_ref(j)^2*original(i).S)
            fprintf('M = %f Nm (Otimizado)\n',archive(i).aero(j,4)*1/2*dat.rho(j)*dat.v_ref(j)^2*archive(i).S*archive(i).mac)
        end
    end
    fprintf('Pontuação: %.6f\n\n\n\n',archive(i).score)
end

%disp(archive(end).aero(1,3)/original.aero(1,3)-1)
disp(archive(end).c_t/archive(end).c_r)
disp(archive(end).b/((archive(end).c_t+archive(end).c_r)/2))


%% Fazer gráfico das médias dos coeficientes [apagar isto depois]
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
%
%saveas(gcf,'aoba.png')


%run_cst_TCC2_3D(v_ex,v_in,dat,[0.5,1],1);
