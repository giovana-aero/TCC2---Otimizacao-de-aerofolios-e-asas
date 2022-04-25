function coo = run_cst_tcc2(dat,v_ex,v_in,op)
    
% Modificação: retirada do fator de escala

% Versão nova da função. Agora admite graus diferentes do polinômio.
% Inputs:
% - Comprimento da corda c
% - Grau do polinômio n
% - Número de pontos np
% - Parâmetros da função shape N1 e N2
% - Vetores com informações do extradorso v_ex
% - Vetores com informações do intradorso v_in
% Dos quais c, n, np, N1 e N2 serão retirados do struct
% Quando op == 1 as coordenadas serão imprimidas

% NOTA: não é necessário imprimir o arquivo de input todas as vezes.
% Pensar nisso uma outra hora

if nargin == 3 % Ignorar a impressão de coordenadas se op sequer estiver no input da função
    op = 0;
end

% pegar informações do struct
c = dat.chord;
n = dat.BPn;
np = dat.np;
N1 = dat.N1;
N2 = dat.N2;
p_op = dat.p_op;

% Inicializar vetores das coordenadas e pesos
if p_op == 1
    x = cosspace(0,c,np);
else    
    x = cosspace_half(0,c,np);
end
y1 = zeros(1,length(x)); y2 = y1;
A1 = zeros(1,n+1); A2 = A1;

% Ler informações dos vetores e calcular pesos
Rle1 = v_ex(1); Rle2 = v_in(1);
beta1 = v_ex(n+1); beta2 = v_in(n+1);
%f1 = v_ex(n+2); f2 = v_in(n+2);
Dz1 = v_ex(n+2); Dz2 = v_in(n+2);
A1(1) = sqrt(2*Rle1/c);
A1(2:n) = v_ex(2:n);
A1(n+1) = tand(beta1) + Dz1/c;
A2(1) = sqrt(2*Rle2/c);
A2(2:n) = v_in(2:n);
A2(n+1) = tand(beta2) + Dz2/c;

%% Extradorso ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for p = 1:np
    
    % Calcular o polinômio de Bernstein
    sum = 0;
    for r = 0:n
        K = factorial(n)/(factorial(r)*factorial(n-r));
        sum = sum + A1(r+1)*K*x(p)^r*(1-x(p))^(n-r);
    end
    
    % Calcular a ordenada com as funções class e shape ao mesmo tempo
    y1(p) = x(p)^N1*(1-x(p))^N2*sum + x(p)*Dz1/c;
end

%% Intradorso ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for p = 1:np
    
    % Calcular o polinômio de Bernstein
    sum = 0;
    for r = 0:n
        K = factorial(n)/(factorial(r)*factorial(n-r));
        sum = sum + A2(r+1)*K*x(p)^r*(1-x(p))^(n-r);
    end
    
    % Calcular a ordenada com as funções class e shape ao mesmo tempo
    y2(p) = -(x(p)^N1*(1-x(p))^N2*sum + x(p)*Dz2/c);
    
end

coo = [flip(x'),flip(y1');
       x(2:end)',y2(2:end)'];
%figure(1),clf
%plot(coo(:,1),coo(:,2)),grid on,axis equal,hold on
%scatter(coo(:,1),coo(:,2))


% Imprimir coordenadas
if op == 1
    ID = fopen('coordenadas.dat','w');
    fprintf(ID,'%f %f\n',coo');
    fclose(ID);
end

 
end
