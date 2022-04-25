function aero = run_xfoil_cst_TCC2(dat,op) 

if op == 1 % Montar o arquivo de input

    % Pegar as informações do código base:
    reynolds = dat.reynolds;
    aoa = dat.aoa;
    iter_sim = dat.iter_sim;
    numNodes = dat.numNodes;
    cases = dat.cases;

    %fclose('all');

    %if (exist('xfoil_input.txt','file'))
    %    delete('xfoil_input.txt');
    %end

    %% Criação do arquivo de input do Xfoil
    fid = fopen('xfoil_input.txt','w');

    % Mudar uma opção dos gráficos do XFOIL (desativar a aparição 
    % da janela com o desenho da simulação)
    fprintf(fid,'PLOP\n');
    fprintf(fid,'G\n\n');

    % Ler coordenadas
    fprintf(fid,'LOAD coordenadas.dat\n\n'); % Nota: retirar a quebra de linha extra 
                                             % caso o perfil tenha um nome
                                             
    % Modificar o número de nós
    if numNodes ~= 0
        fprintf(fid,'PPAR\n');
        fprintf(fid,'N %s\n', num2str(numNodes));
        fprintf(fid,'\n\n');
    end
                                            
    % Simulação
    fprintf(fid,'OPER\n');
    if cases == 1 % Apenas uma condição de voo
        
        fprintf(fid,'VISC %d\n', reynolds(1)); % Aplicar o modo viscoso e determinar número de Reynolds 
        fprintf(fid,'PACC\n');                 % Estabelecer arquivo de output
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
            if i == 1 || i~= 1 && iter_sim(i) ~= iter_sim(i-1)
                fprintf(fid,'ITER %d\n', iter_sim(i)); % Mudar o número de iterações    
            end
            fprintf(fid,'ALFA %f\n', aoa(i));        % Estabelecer ângulo de ataque
            fprintf(fid,'PACC\n'); % Fechar a polar para começar a próxima simulação
            
        end
    end

    % Fechar arquivo
    fprintf(fid,'\nQUIT\n');
    fclose(fid);

    
else % Simular os aerofólios
    cases = dat.cases;
    
    % Apagar arquivos caso existam
    for i = 1:cases
        a_polar = ['polar' num2str(i) '.txt'];
        if (exist(a_polar,'file'))
            delete(a_polar);
        end
    end

    % Executar XFOIL com o arquivo de entrada
    cmd = 'xfoil.exe < xfoil_input.txt';
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
            
            aero(i,:) = [CL,CD,CL/CD,CM];
            
        else
            aero = 'n'; break % Ignorar o aerofólio se ao menos uma das simulações não convergir
        end

        
    end


end