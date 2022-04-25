%% Algoritmo genético - Perfis NACA 4 dígitos 
% Versão TCC 2

% Modificações a serem feitas:
% - Na checagem de erros, pôr código que troca o número de iterações do XFOIL pro número padrão quando 
%   valores nulos forem inseridos



%%

clear,clc
fclose('all');
tic 
[y2,Fs2] = audioread('C:\Users\Guga Weffort\Documents\MATLAB\504 Finding a Treasure Box.wav');


% Parâmetros do algoritmo
dat.N = 500;                              % Número de indivíduos na população
dat.mu = 0.05;                           % Probabilidade de mutação
dat.iter = 5;                           % Número de iterações
dat.elite = 1;                         % Aplicar elitismo?
dat.subs = 1;                          % Substituir aerofólios sem resultados? 

% Parâmetros da geometria
dat.m_ext = [0,9];
dat.p_ext = [1,9];
dat.t_ext = [10,30];

% Parâmetros das simulações
dat.cases = 2;                          % Número de condições de voo a serem analisadas
dat.reynolds = [1e6,1e6,1e6];           % Valores dos números de Reynolds para as simulações
dat.aoa = [0,4,0];                     % Ângulos de ataque
dat.iter_sim = [10,10,10];             % Números de iterações no XFOIL
dat.coeff_op = ['o','^','!','c';       % Uma linha para cada condição de voo
                '!','!','!','!';
                '!','!','!','!'];
dat.coeff_val = [0.2,0.021,90,0;
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
% A matriz dat.coeff_F dá os pesos de cada função objetiva

% Checagem de erros
dat = error_check_naca4_TCC2(dat);

% Gerar população inicial
empty.m = [];
empty.p = [];
empty.t = [];
empty.aero = []; % Terá o mesmo formato que a matriz coeff_op
empty.score = 0;

pop = repmat(empty,dat.N,1);
chi = pop;
    
% Gerar numeração aleatória
disp('<< Geração da população inicial >>')
for i = 1:dat.N
    disp(['Indivíduo ' num2str(i)])
    pop(i).m = randi(dat.m_ext); % Curvatura máxima
    if pop(i).m == 0 % Perfis simétricos
        pop(i).p = 0;
    else % Perfis assimétricos
        pop(i).p = randi(dat.p_ext);
    end
    pop(i).t = randi(dat.t_ext); % Espessura máxima
end

% Gerar struct que guarda o melhor perfil de cada geração
archive = repmat(empty,dat.iter,1);

%% Loop principal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for loop = 1:dat.iter

    % Simular os perfis e obter dados
    select = ones(1,dat.N);
    disp('<< Simulação dos aerofólios >>')
    for i = 1:dat.N
        disp(['Indivíduo ' num2str(i)])
        naca_num = strcat(num2str(pop(i).m),num2str(pop(i).p),num2str(pop(i).t));
        pop(i).aero = run_xfoil_naca4_TCC2(naca_num,dat);clc % Simulação
        
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
    pop = fitness_naca4(pop,dat,select2);
    
    % Isto serve pra pôr todas as pontuações em um vetor
    weights = [pop.score];
    
    % Guardar o melhor perfil de cada iteração
    [~,pos] = max(weights);
    archive(loop) = pop(pos);
    
    % Mostrar o melhor perfil
    figure(1),clf
    plot_airfoil_naca4_TCC2(pop(pos),loop)
    
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
        if pop(par(1)).m == 0 || pop(par(2)).m == 0 % Caso um deles seja simétrico
            n = 1;
        else
            n = randi(3);
        end
        
        switch n
            case 1
                chi(c).m = pop(par(1)).m;
                chi(c).p = pop(par(1)).p;
                chi(c).t = pop(par(2)).t;
                chi(c+1).m = pop(par(2)).m;
                chi(c+1).p = pop(par(2)).p;
                chi(c+1).t = pop(par(1)).t;
                
            case 2
                chi(c).m = pop(par(1)).m;
                chi(c).p = pop(par(2)).p;
                chi(c).t = pop(par(1)).t;
                chi(c+1).m = pop(par(2)).m;
                chi(c+1).p = pop(par(1)).p;
                chi(c+1).t = pop(par(2)).t;
                
            case 3
                chi(c).m = pop(par(2)).m;
                chi(c).p = pop(par(1)).p;
                chi(c).t = pop(par(1)).t;
                chi(c+1).m = pop(par(1)).m;
                chi(c+1).p = pop(par(2)).p;
                chi(c+1).t = pop(par(2)).t;
                
        end
        if (chi(c).m == 0 && chi(c).p ~= 0) || (chi(c).m ~= 0 && chi(c).p == 0) || (chi(c+1).m == 0 && chi(c+1).p ~= 0) || (chi(c+1).m ~= 0 && chi(c+1).p == 0)
            error('problema na configuração do perfil simétrico')
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
            
            
            % Código da mutação referente apenas aos perfis NACA 4 dígitos
            % começa aqui
            if chi(s).m == 0 % Perfis simétricos
                n = randi(2);
                switch n
                    case 1 % Transformar em um perfil assimétrico
                        chi(s).m = randi([1,dat.m_ext(2)]);
                        chi(s).p = randi(dat.p_ext);
                    case 2 % Mudar apenas a espessura
                        chi(s).t = randi(dat.t_ext);
                end
            else % Perfis assimétricos
                n = randi(4);
                switch n
                    case 1 % Mudar apenas a curvatura máxima
                        chi(s).m = randi([1,dat.m_ext(2)]);
                    case 2 % Mudar apenas o local da curvatura máxima
                        chi(s).p = randi(dat.p_ext);                        
                    case 3 % Mudar apenas a espessura
                        chi(s).t = randi(dat.t_ext);
                    case 4 % Transformar em um perfil simétrico
                        chi(s).m = 0;
                        chi(s).p = 0;
                end
            end
                
                
            % Código da mutação referente apenas aos perfis NACA 4 dígitos
            % termina aqui
            
            if chi(s).m == 0 && chi(s).p ~= 0 || chi(s).m ~= 0 && chi(s).p == 0 
                error('problema na configuração do perfil simétrico')
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



%% Final ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Mostrar qual foi o perfil resultante, as condições da simulação e o
% número de iterações. Também traçar um gráfico do perfil

clc,%sound(y2,Fs2)
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
    
    % Trocar separador decimal
    xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
    new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
    set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

end

%% Pegar o struct de arquivo e imprimir todos
for i = 1:length(archive)
    disp(['<< Iteração ' num2str(i) ' >>'])
    fprintf('NACA %d%d%d\n',archive(i).m,archive(i).p,archive(i).t)
	disp('- Dados aerodinâmicos -')
    for j = 1:dat.cases
        fprintf('Condição de voo %d: ',j),disp(archive(i).aero(j,:))
    end
    fprintf('Pontuação: %.6f\n\n\n',archive(i).score)%,disp(archive(i).score)
    
end



