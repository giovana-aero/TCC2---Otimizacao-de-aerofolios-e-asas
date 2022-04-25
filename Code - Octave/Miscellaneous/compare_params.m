clc,clear

%v1 = [7.1939e+01,7.5612e+01,7.6268e+01,8.1434e+01,8.1434e+01];
%v2 = [7.1844e+01,7.2615e+01,7.9603e+01,8.3381e+01,8.4203e+01];
%v3 = [8.0220e+01,8.0220e+01,8.6911e+01,8.6911e+01,8.6911e+01];
v1 = [7.3524e+01,7.1655e+01,6.8510e+01,7.4659e+01,7.7161e+01,7.7718e+01,8.3665e+01,8.4542e+01,8.7220e+01,8.2440e+01];
v2 = [6.8139e+01,7.0118e+01,8.3462e+01,8.4707e+01,8.5888e+01,8.9817e+01,8.9817e+01,8.9817e+01,9.0914e+01,9.3408e+01];


iter = 1:10;

figure(1),clf
plot(iter,v1,'k-*'),grid on,hold on
plot(iter,v2,'r-*')
%plot(iter,v3,'b-*')
xlabel('Iteração'),ylabel('L/D')
%legend('Mutação 5%','Mutação 10%','Mutação 20%')
legend('Elitismo desativado','Elitismo ativado',"location","southeast")
saveas(gcf,'aoba.png')