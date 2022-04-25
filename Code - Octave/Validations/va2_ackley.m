%% Validação 2 - Função de Ackley


clc,clear,clf

% Parâmetros do algoritmo
dat.n = 100;                            % Número de indivíduos na população (deve ser par por hora)
dat.iter = 200;                         % Número de iterações
dat.t = 20;                               % Dimensão da função
dat.mu = 0.05;                           % Probabilidade de mutação
dat.mu_number = 1/dat.mu;                % Vetor usado na fase de mutação
dat.best = zeros(1,dat.iter);            % Vetor usado no plot final (contém o valor de função que mais se aproxima do mínimo)
dat.convergence_op = 1;

%% Gerar população inicial

% Template
empty.gen = [];  % Genótipo
empty.score = 0; % Pontuação (quanto maior, melhor)

% População inicial
pop = repmat(empty,dat.n,1);
children = pop;

% Atribuir caracteres aleatórios aos genes de cada indivíduo
for i = 1:dat.n
    for j = 1:dat.t
        pop(i).gen(j) = rand*32.768*randi([-1 1]);
    end
end

% condition = 0; 
convergence = 0;
%% loop principal
for loop = 1:dat.iter
    
    % Avaliar cada indivíduo de acordo com a função e atribuir pontuação
    for i = 1:dat.n
        
        sum1 = 0;
        sum2 = 0;
        for j = 1:dat.t
            sum1 = pop(i).gen(j)^2. + sum1;
            sum2 = cos(2*pi*pop(i).gen(j)) + sum2;
        end
        
        f_tgt = -20*exp(-1/5*sqrt(dat.t^-1*sum1)) - exp(dat.t^-1*sum2) + 20 + exp(1);
        
        if f_tgt == 0
            condition = 1;
            break
            
        else
            pop(i).score = f_tgt^-1;
            
        end
    end
    
    
    % Isto serve pra pôr todas as pontuações em um vetor
    weights = zeros(1,dat.n);
    for i = 1:dat.n
        weights(i) = pop(i).score;
    end
    
    % Pegar o melhor de cada iteração
    [~,pos] = max(weights);
    sum1 = 0;
    sum2 = 0;
    for j = 1:dat.t
        sum1 = pop(pos).gen(j)^2. + sum1;
        sum2 = cos(2*pi*pop(pos).gen(j)) + sum2;
    end
    f_tgt = -20*exp(-1/5*sqrt(dat.t^-1*sum1)) - exp(dat.t^-1*sum2) + 20 + exp(1);
    dat.best(loop) = f_tgt;
    fprintf('Iteraçao %d - Indivíduo %d, ',loop,pos)
    fprintf('com f = %f\n',f_tgt)
    
    
    % Checar convergência
    if dat.convergence_op == 1
        check = (weights == pop(1).score);
        sum = 0;
        for i = 1:size(check,2)
            if check(i) == 1
                sum = sum + 1;
                convergence = 1;
            end
            
        end
        
        if sum>0.9*dat.n
            break
        end
        
    end
    
    
    % Escolher membros da população pra reprodução
    c=1;
    for f = 1:dat.n/2
        
        % Isto seleciona dois pais por meio de uma seleção via roleta
        % (indivíduos com pesos maiores têm  mais chance de serem selecionados)
        par = [0 0]; % Vetor que indica a numeração dos pais escolhidos
        for u = 1:2
            accumulation = cumsum(weights);
            p = rand() * accumulation(end);
            chosen_index = -1;
            for index = 1 : length(accumulation)
                if (accumulation(index) > p)
                    chosen_index = index;
                    break;
                end
            end
            %choice = chosen_index;
            par(1,u) = chosen_index;
        end
        
        
        % Crossover, gerando dois filhos de cada dois pais
        corte = randi(dat.t);
        children(c).gen = [pop(par(1,1)).gen(1:corte) pop(par(1,2)).gen(1+corte:end)];
        children(c+1).gen = [pop(par(1,2)).gen(1:corte) pop(par(1,1)).gen(1+corte:end)];
        c = c + 2;
        
        
    end
    
    
    % Mutação
    for i = i:dat.n
        
        if randi([1 dat.mu_number]) == 1
            
            start = randi(dat.t);
            range = randi([start dat.t]);
            array = zeros(1,range-start+1);
            
            k = 1;
            for j = start:range
                children(i).gen(j) = rand*32.768*randi([-1 1]);
                k = k + 1;
            end
            
        end
    end
    
    
    
    % Substituir a população inicial pelos filhos
    pop = children;
    
end



if convergence == 1
    disp('Solução convergida')
end

plot(1:loop,dat.best(1:loop)),grid on
xlabel('Iteração'),ylabel('f(x)')


% sum = 0;
% for j = 1:dat.t
%     sum = pop(best2).gen(j)^2. + sum;
% end
% f_tgt = sum;
% fprintf('Valor da função: %f\n',f_tgt)



