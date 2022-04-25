function [v_ex,v_in] = reverse_CST_LS(coo,n)
% Função: transformar um conjunto de coordenadas de aerofólio em vetores CST via
% método dos mínimos quadrados

% Definir a função dos coeficientes binomiais 
function K = binom(r,n)
    K = factorial(n)/(factorial(r)*factorial(n-r));
end

% Transladar as coordenadas caso a menor ordenada não seja igual a zero
[val,pos] = min(coo(:,1));
if val ~= 0
    coo(:,1) = coo(:,1) - val;
    coo(:,2) = coo(:,2) - coo(pos,2);
end
if coo(pos,2) ~= 0
    coo(:,2) = coo(:,2) - coo(pos,2);
end

% Tomar o ponto de separação das superfícies no bordo de ataque como a ordenada 
% mínima para facilitar o funcionamento do CST
ex = flip(coo(1:pos,:)); 
in = coo(pos:end,:);
c = coo(1,1);
np = 80;

% Encontrar as distâncias verticais do bordo de fuga
Dz1 = (ex(end,2) - coo(pos,2));
Dz2 = -(in(end,2) + coo(pos,2));

% Montar o sistema de equações(extradorso) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
M = zeros(n+2);
R = zeros(n+2,1);

% Pegar as ordenadas e abscissas do extradorso para facilitar
ex_x = ex(:,1);
ex_y = ex(:,2);

for i = 1:n+1

    K_a = binom(i-1,n);
    g_a = (ex_x).^(0.5).*(1-ex_x).*K_a.*(ex_x).^(i-1).*(1-ex_x).^(n-(i-1));
        
    for j = 1:n+1
        K_b = binom(j-1,n);
        g_b = (ex_x).^(0.5).*(1-ex_x).*K_b.*(ex_x).^(j-1).*(1-ex_x).^(n-(j-1));
        
        M(i,j) = sum(g_a.*g_b);
    end
    
    % Termo com o bordo de fuga
    g_b = ex_x./c;
    M(i,j+1) = sum(g_a.*g_b);
    
    R(i,1) = sum(ex_y.*g_a);
end

% Últimas linhas das matrizes M e R
g_a = ex_x./c;
for j = 1:n+1
    K_b = binom(j-1,n);
    g_b = (ex_x).^(0.5).*(1-ex_x).*K_b.*(ex_x).^(j-1).*(1-ex_x).^(n-(j-1));
    
    M(end,j) = sum(g_a.*g_b);
end

% Termo com o bordo de fuga
g_b = ex_x./c;
M(end,j+1) = sum(g_a.*g_b);
R(end,1) = sum(ex_y.*g_a);

% Solucionar o sistema linear
A1 = linsolve(M,R)';

% Montar o sistema de equações (intradorso) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
M = zeros(n+2);
R = zeros(n+2,1);

% Pegar as ordenadas e abscissas do extradorso para facilitar
in_x = in(:,1);
in_y = -in(:,2); % Trocar o sinal aqui pra facilitar as contas

for i = 1:n+1

    K_a = binom(i-1,n);
    g_a = (in_x).^(0.5).*(1-in_x).*K_a.*(in_x).^(i-1).*(1-in_x).^(n-(i-1));
        
    for j = 1:n+1
        K_b = binom(j-1,n);
        g_b = (in_x).^(0.5).*(1-in_x).*K_b.*(in_x).^(j-1).*(1-in_x).^(n-(j-1));
        
        M(i,j) = sum(g_a.*g_b);
    end
    
    % Termo com o bordo de fuga
    g_b = in_x./c;
    M(i,j+1) = sum(g_a.*g_b);
    
    R(i,1) = sum(in_y.*g_a);
end

% Últimas linhas das matrizes M e R
g_a = in_x./c;
for j = 1:n+1
    K_b = binom(j-1,n);
    g_b = (in_x).^(0.5).*(1-in_x).*K_b.*(in_x).^(j-1).*(1-in_x).^(n-(j-1));
    
    M(end,j) = sum(g_a.*g_b);
end

% Termo com o bordo de fuga
g_b = in_x./c;
M(end,j+1) = sum(g_a.*g_b);
R(end,1) = sum(in_y.*g_a);

% Solucionar sistema linear
A2 = linsolve(M,R)';

%% Montar vetores
v_ex = [A1(1)^2*c/2,A1(2:n),atand(A1(n+1)-Dz1/c),Dz1];
v_in = [A2(1)^2*c/2,A2(2:n),atand(A2(n+1)-Dz2/c),Dz2];


end