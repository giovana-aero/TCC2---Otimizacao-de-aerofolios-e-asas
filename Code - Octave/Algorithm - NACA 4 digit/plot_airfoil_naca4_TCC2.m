function plot_airfoil_naca4_v2(pop,loop)

Chord = 1;
x = 0:0.01:Chord;    
m = pop.m/100;
p = pop.p*Chord/10;
t = pop.t/100;
%aero = pop.aero;

if m == 0
    Symm = 1;
else
    Symm = 0;
end

if Symm == 1
    y_upper = 5*t*Chord*(0.2969*sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord).^2+0.2843*(x/Chord).^3-0.1015*(x/Chord).^4); 
    y_lower = -y_upper;
    x_upper = x;
    x_lower = x;
else
    y_camber = zeros(1,length(x));
    dy_camber = zeros(1,length(x));
    for i = 1:length(x)
        if x(i)/Chord<=p
            y_camber(i) = m*x(i)/p^2*(2*p-x(i)/Chord);
            dy_camber(i) = 2*m/p^2*(p-x(i)/Chord);
        else
            y_camber(i) = m*(Chord-x(i))/(1-p)^2*(1+x(i)/Chord-2*p);
            dy_camber(i) = 2*m/(1-p)^2*(p-x(i)/Chord);
        end
    end
    y_t = 5*t*Chord*(0.2969*sqrt(x/Chord)-0.126*(x/Chord)-0.3516*(x/Chord).^2+0.2843*(x/Chord).^3-0.1015*(x/Chord).^4); 
    theta = atan(dy_camber);
    x_upper = x-y_t.*sin(theta);
    x_lower = x+y_t.*sin(theta);
    y_upper = y_camber+y_t.*cos(theta);
    y_lower = y_camber-y_t.*cos(theta);
end

% Traçar o perfil
plot([flip(x_upper) x_lower],[flip(y_upper) y_lower])
axis equal,grid on
title(['Iteração ' num2str(loop) ': NACA ' num2str(pop.m) num2str(pop.p) num2str(pop.t)])

xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)


end