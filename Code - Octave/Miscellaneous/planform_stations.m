clc,clear

b = 10;
cr = 3;
ct = 1;

n = 5;


% Primeira equação
B1 = cr/2;
A1 = (ct/2-B1)*2/b;

% Segunda equação


x2 = linspace(0,1,n+2)*b/2;
x1 = -flip(x2);
%cm = flip(2*x1*cr/b);
y = flip(A1*x2+B1);

% Ponta esquerda
figure(1),clf
%plot([x1(1),ct/2],[x1(1),-ct/2],'k')
plot([x1(1),x1(1)],[ct/2,-ct/2],'k','linewidth',3),hold on,axis equal
% Meio esquerdo
for i = 2:n+1
    plot([x1(i),x1(i)],[y(i),-y(i)],'r','linewidth',3)
end
% Raiz
plot([0,0],[-cr/2,cr/2],'k','linewidth',3)
% Meio direito
y = flip(y);
for i = 2:n+1
    plot([x2(i),x2(i)],[y(i),-y(i)],'r','linewidth',3)
end
% Ponta direita
plot([x2(end),x2(end)],[ct/2,-ct/2],'k','linewidth',3)

%saveas(gcf,'aoba.png')