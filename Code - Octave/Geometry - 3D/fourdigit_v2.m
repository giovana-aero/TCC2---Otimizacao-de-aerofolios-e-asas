function coo = fourdigit_v2(naca_num,np,x_op)
    
% Entradas:
% Nome do perfil (vetor 1x3: [m,p,t])
% np - número de pontos (abscissa)
% x_op - opção de geração de pontos (1 pra cosspace, outro valor pra cosspace_half)

%if nargin == 4,x_op=1;end

if nargin == 2 || x_op == 1  % Opção padrão
    x = cosspace(0,1,np);
else
    x = cosspace_half(0,1,np);
end

m = naca_num(1)/100;
p = naca_num(2)/10;
t = naca_num(3)/100;

% Distribuição de espessura
% (o último valor muda de -0.1015 para -0.1036 no caso de bordo de fuga
% fechado)
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1036;
yt = 5*t*(a0*x.^0.5+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

if m == 0 && p == 0 % Perfil simétrico
    xU = x;
    yU = yt;
    xL = x;
    yL = -yt;
    
else  % Perfil com curvatura
    
    % Distribuição de curvatura
    yc = zeros(1,size(x,2));
    for i = 1:size(x,2)
        if x(i) >= 0 && x(i) <= p
            yc(i) = m/p^2*(2*p*x(i)-x(i)^2);
        elseif x(i) >= p && x(i) <= 1
            yc(i) = (m/(1-p)^2)*((1-2*p)+2*p*x(i)-x(i)^2);
        end
    end
    
    derivada = zeros(1,size(x,2));
    for i = 1:size(x,2)
        if x(i) >= 0 && x(i) <= p
            derivada(i) = 2*m/p^2*(p-x(i));
        elseif x(i) >= p && x(i) <= 1
            derivada(i) = 2*m/(1-p)^2*(p-x(i));
        end
    end
    
    theta = atan(derivada);
    xU = x - yt.*sin(theta);
    yU = yc + yt.*cos(theta);
    xL = x + yt.*sin(theta);
    yL = yc - yt.*cos(theta);
    
end

% Coordenadas
coo = [flip(xU'),flip(yU');
xL(2:end)',yL(2:end)'];
       
end