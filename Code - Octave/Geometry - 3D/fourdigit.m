function [xU,yU,xL,yL] = fourdigit(x,m,p,t)

m = m/100;
p = p/10;
t = t/100;

% Simétrico

% distribuição de espessura
% (o último valor muda de -0.1015 para -0.1036 no caso de bordo de fuga
% fechado)
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1036;

yt = 5*t*(a0*x.^0.5+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);


if m == 0 && p == 0    % perfil simétrico
    xU = x;
    yU = yt;
    xL = x;
    yL = -yt;
    
    
else  % perfil com curvatura
    
    % distribuição de curvatura
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
    
    
    % coordenadas
    xU = x - yt.*sin(theta);
    yU = yc + yt.*cos(theta);
    xL = x + yt.*sin(theta);
    yL = yc - yt.*cos(theta);
    
    
    
end

%plot(xU,yU,xL,yL)
%axis square, set(gca,'ylim',[-0.5 0.5],'xlim',[0 1])
%axis equal,grid on



end