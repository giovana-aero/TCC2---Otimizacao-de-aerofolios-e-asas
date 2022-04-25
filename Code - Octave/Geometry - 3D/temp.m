clc,clear

c_r = 2;
c_t = 1;
op = 0:0.2:1;

%coo_r = 'coordenadas.dat';
%coo_t = 'coordenadas_2.dat';
%coo_r = dlmread(coo_r)*c_r;
%coo_t = dlmread(coo_t)*c_t;
x = cosspace_half(0,1,80);
[xU,yU,xL,yL] = fourdigit(x,9,4,16); coo_r = [flip(xU'),flip(yU');xL(2:end)',yL(2:end)']*c_r;
[xU,yU,xL,yL] = fourdigit(x,2,4,12); coo_t = [flip(xU'),flip(yU');xL(2:end)',yL(2:end)']*c_t;



% Traçar gráficos dos perfis originais
figure(1),clf
plot(coo_r(:,1),coo_r(:,2)),axis equal,hold on,grid on
plot(coo_t(:,1),coo_t(:,2))


for i = 1:length(op)
    coo = airfoil_interpolation(coo_r,coo_t,op(i));
    plot(coo(:,1),coo(:,2),'--')
end
