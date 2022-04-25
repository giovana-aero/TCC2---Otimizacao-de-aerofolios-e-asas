function run_CST_old(c,x_c,v_ex,v_in)

n = 80;
N = 4;   % Grau do polinômio

Rle = v_ex(1);

% Pesos da geometria - extradorso
A1  = zeros(1,N+1);
A1(2) = v_ex(2);
A1(3) = v_ex(3);
A1(4) = v_ex(4);
beta1 = v_ex(5);
f1 = v_ex(6);
delta_z1 = v_ex(7);

% Pesos da geometria - intradorso
A2 = zeros(1,N+1);
A2(2) = v_in(2);
A2(3) = v_in(3);
A2(4) = v_in(4);
beta2 = v_in(5);
f2 = v_in(6);
delta_z2 = v_in(7);

%% Funções (extradorso)

% Class
C1 = x_c.^0.5.*(1-x_c);

% Shape
S1 = zeros(1,n);
A1(1) = sqrt(2*Rle/c);
A1(end) = tand(beta1) + delta_z1/c;
S1(1) = A1(1);
S1(end) = A1(end);

for i = 1:n
    
    soma = zeros(1,length(A1));
    
    for j = 0:N
        K = factorial(N)/(factorial(j)*factorial(N-j));
        soma(j+1) = A1(j+1)*K*x_c(i)^j*(1-x_c(i))^(N-j);
    end
    
    S1(i) = sum(soma);
    
end

z_c1 = C1.*S1 + x_c.*delta_z1;

%% Funções (intradorso)

% Class
C2 = x_c.^0.5.*(1-x_c);

% Shape
S2 = zeros(1,n);
A2(1) = sqrt(2*Rle/c);
A2(end) = tand(beta2) + delta_z2/c;
S2(1) = A2(1);
S2(end) = A2(end);

for i = 1:n
    
    soma = zeros(1,length(A2));
    
    for j = 0:N
        K = factorial(N)/(factorial(j)*factorial(N-j));
        soma(j+1) = A2(j+1)*K*x_c(i)^j*(1-x_c(i))^(N-j);
    end
    
    S2(i) = sum(soma);
    
end

z_c2 = C2.*S2 + x_c.*delta_z2;

% plot(x_c,z_c1*f1,x_c,-z_c2*f2),axis equal

%% Gerar arquivo de coordenadas

if (exist('coordenadas.dat','file'))
    delete('coordenadas.dat');
end

file = fopen('coordenadas.dat','w');
%fprintf(file,'Perfil CST teste\n');
for i = 0:n-1
    fprintf(file,'%f %f\n',x_c(end-i),z_c1(end-i)*f1);
end
for i = 2:n
    fprintf(file,'%f %f\n',x_c(i),-z_c2(i)*f2);
end

fclose(file);

end