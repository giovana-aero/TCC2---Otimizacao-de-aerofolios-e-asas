clc,clear

nome = 'wing_test.tri';

coo = dlmread(nome);
figure(1),clf
scatter3(coo(:,1),coo(:,2),coo(:,3)),grid on,axis equal

p1 = [1,1];
p2 = [2,2];
plot3([1,2,0],[3,4,0])
xlabel('x'),ylabel('y'),zlabel('z')