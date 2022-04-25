function check = quality(coo,dat)
% Checagem de qualidade dos perfis CST

np = dat.np;
c = dat.chord;

% Separar as superfícies
ex = flip(coo(1:np,2)); in = coo(np:(np*2)-1,2);
% Definir a espessura
thickness = ex - in;

check = 1;
while 1
    
    % 1 - Checar se há interseção
    if sum(thickness(2:end-1)<=0) > 0
        check = 0;disp('erro1'),break
    end
    
    % 2 - Checar as inclinações do extradorso
    % A inclinação do extradorso não pode mudar duas ou mais vezes
    slope = zeros(1,length(ex));
    for i = 1:(length(ex)-1)
        slope(i) = ex(i+1) - ex(i);
    end
    v = slope > 0;
    
    % Isto detecta quantas mudanças de inclinação existem
    counter = 0; num = 0;
    for i = 1:length(v)
        if v(i) == num
            counter = counter + 1;
            if mod(counter,2) == 0
                num = 0;
            else
                num = 1;
            end
        end
    end
    
    if counter >= 2
        check = 0; disp('erro2'),break
    end
    
    % 3 - Checar as inclinações do intradorso
    % A inclinação do intradorso não pode mudar três ou mais vezes
    slope = zeros(1,length(in));
    for i = 1:(length(in)-1)
        slope(i) = in(i+1) - in(i);
    end
    v = slope < 0;
    
    counter = 0; num = 0;
    for i = 1:length(v)-1
        if v(i) == num
            counter = counter + 1;
            if mod(counter,2) == 0
                num = 0;
            else
                num = 1;
            end
        end
    end
    
    if counter >= 3
        check = 0; disp('erro3'),break
    end
    
    % 4 - Garantir que o perfil não seja fino demais
    if mean(thickness)/c < 0.04
        check = 0; disp('erro4'),break
    end
    
	%% 5 - Garantir que o bordo de fuga não seja pontiagudo demais
    %if mean(ex(end-5:end-1) - in(end-5:end-1))/c < 7e-3
    %    check = 0; break
    %end
	
    break
end