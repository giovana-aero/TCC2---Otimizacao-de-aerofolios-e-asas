clc,clear

op = 0.5;
naca1 = [4,4,12];
naca2 = [0,0,10];
naca3 = [naca1(1)*(1-op)+naca2(1)*op,naca1(2)*(1-op)+naca2(2)*op,naca1(3)*(1-op)+naca2(3)*op];
n = 100;


%naca1 = str(naca1[0]) + str(naca1[1]) + str(naca1[-1])
%naca2 = str(naca2[0]) + str(naca2[1]) + str(naca2[-1])
%naca3 = str(naca3[0]) + str(naca3[1]) + str(naca3[-1])

%x1,y1 = naca4(naca1,n)
%x2,y2 = naca4(naca2,n)
%x3,y3 = naca4(naca3,n)
%
%coo1 = np.hstack((np.array([x1]).transpose(),np.array([y1]).transpose()))
%coo2 = np.hstack((np.array([x2]).transpose(),np.array([y2]).transpose()))
%coo3 = np.hstack((np.array([x3]).transpose(),np.array([y3]).transpose()))

coo1 = fourdigit_v2(naca1,100);
coo2 = fourdigit_v2(naca2,100);
coo3 = fourdigit_v2(naca3,100);

coo_intp = airfoil_interpolation(coo1,coo2,op,0,0);

figure(1),clf
plot(coo1(:,1),coo1(:,2),'k',coo2(:,1),coo2(:,2),'k--',coo_intp(:,1),coo_intp(:,3),'r'),hold on,grid on
%plot(coo3(:,1),coo3(:,2),'--')'
axis equal
%legend('raiz','ponta','interpolado (função)','interpolado (manipulando os números)')
legend('Raiz','Ponta','Interpolado')

set(gca,'ylim',[-0.2,0.3])
set(gca,'xlim',[0,1])

% saveas(gcf,'aoba.png')
