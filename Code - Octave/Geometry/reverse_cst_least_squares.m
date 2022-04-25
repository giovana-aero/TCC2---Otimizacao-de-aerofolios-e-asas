%% Este código lê as coordenadas de algum aerofólio e fornece os parâmetros 
%% CST que geram a mesma geometria por meio do método dos mínimos quadrados

% As coordenadas de entrada DEVEM estar configuradas como:
% Bordo de fuga -> bordo de ataque -> bordo de fuga


clc,clear

% Definir a função dos coeficientes binomiais 
function K = binom(r,n)
    K = factorial(n)/(factorial(r)*factorial(n-r));
end

% Qual o grau do polinômio de Bernstein?
% Número de variáveis da função shape (ou seja, não considera o delta z) é igual
% ao grau do polinômio mais um
n = 6;

% Precisão dos valores nos vetores (número de casas decimais)
pre = 6;

% Configuração pra geração de pontos
p_op = 1;

% Ler coordenadas de arquivos de texto
coo = dlmread('coordenadas.dat');
%coo = dlmread('naca2412.txt');
%coo = dlmread('naca5311.txt');
%coo = dlmread('naca0012.txt');
%coo = dlmread('naca0010.txt');
%coo = dlmread('s1223.txt');
%coo = dlmread('fx61-180.txt');
%coo = dlmread('eppler544.txt');
%coo = dlmread('nplx.txt');
%coo = dlmread('raf15.txt');
%coo = dlmread('saratov.txt');
%coo = dlmread('whitcomb.txt');
%coo = dlmread('naca0006.txt');
%coo = dlmread('fx 61-163.txt');
%coo = dlmread('fx 60-126.txt');
%coo = dlmread('naca650006.txt');
%coo = dlmread('NACA 65A209.txt');
%coo = dlmread('NACA 65A210.5.txt');
%coo = dlmread('NACA 65A212.txt');
%coo = dlmread('rae2822.txt');
%coo = dlmread('sd7062.txt');
%coo = dlmread('naca23012.txt');
%coo = dlmread('naca4412.txt');
%%coo = dlmread('naca2312.txt');
%coo = dlmread('naca2411.txt');
coo = converter_function_xfoil(coo);

% Obter coordenadas CST
%v_ex = [0.0200, 0.2700, 0.2400, 0.0700, 21.0000, 0.0000];
%v_in = [0.0200, -0.0000, 0.0200, 0.0500, 0.0000, 0.0000];
%dat.chord = 1;
%dat.BPn = length(v_ex)-2;
%dat.np = 80;
%dat.N1 = 0.5;
%dat.N2 = 1;
%dat.p_op = p_op;
%coo = run_cst_tcc2(dat,v_ex,v_in,1);

% Mostrar perfil original e sua distribuição de pontos
%figure(1),clf
%plot(coo(:,1),coo(:,2)),hold on,grid on,axis equal
%scatter(coo(:,1),coo(:,2)),pause(3)

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


% Pós-processamento ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Imprimir as informações
v_ex = [A1(1)^2*c/2,A1(2:n),atand(A1(n+1)-Dz1/c),Dz1];
v_in = [A2(1)^2*c/2,A2(2:n),atand(A2(n+1)-Dz2/c),Dz2];

% Extradorso
fprintf(['v_ex = [%.' num2str(pre) 'f, '], v_ex(1))
for j = 2:(length(v_ex)-2)
    fprintf(['%.' num2str(pre) 'f, '],v_ex(j))
end
fprintf(['%.' num2str(pre) 'f, '], v_ex(end-1))
fprintf(['%.' num2str(pre) 'f];\n'], v_ex(end))

% Intradorso
fprintf(['v_in = [%.' num2str(pre) 'f, '], v_in(1))
for j = 2:(length(v_in)-2)
    fprintf(['%.' num2str(pre) 'f, '],v_in(j))
end
fprintf(['%.' num2str(pre) 'f, '], v_in(end-1))
fprintf(['%.' num2str(pre) 'f];\n'], v_in(end))

  
% Comparar aerofólios
dat.chord = c;
dat.BPn = n;
dat.np = np;
dat.N1 = 0.5;
dat.N2 = 1;
dat.p_op = p_op;
coo2 = run_cst_tcc2(dat,v_ex,v_in);
figure(1),clf
scatter(coo2(:,1),coo2(:,2)),grid on,hold on,axis equal % Perfil CST calculado
plot(coo(:,1),coo(:,2),'k') % Perfil original
plot(coo2(:,1),coo2(:,2),'r--')
legend('Calculado','Original','Calculado (linha)')

% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)


%clc
%x = log((2*factorial(5)+factorial(8-1))^(sqrt(9))+factorial(4)+factorial(factorial(3)))/sqrt(67);
%fprintf('%.60f\n%.60f\n',pi,x)
