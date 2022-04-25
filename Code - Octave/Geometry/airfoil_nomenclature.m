% Script utilizado pra fazer uma figura pro documento do TCC


clc,clear

%name = 'perfil_cst.txt';
name = 'naca2412.txt';
coo = dlmread(name);
num = 41;

ex = flip(coo(1:num,2));
in = coo(num:end,2);
%camber = (z_c1-z_c2)/2;
camber = (ex + in)/2;

n = 1;
figure(1),clf
plot(coo(:,1),coo(:,2),'k','linewidth',n),hold on,axis equal
%scatter(coo(:,1),coo(:,2),'k')
plot([0,1],[0,0],'b','linewidth',n)
plot(flip(coo(1:num,1)),camber,'r--','linewidth',n)
%plot()


