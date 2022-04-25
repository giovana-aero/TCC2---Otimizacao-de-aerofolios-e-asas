clc,clear

naca1 = [1,1,11];
naca2 = [3,3,11];
%naca3 = [2,3,11];

naca1_str = [num2str(naca1(1)),num2str(naca1(2)),num2str(naca1(3))];
naca2_str = [num2str(naca2(1)),num2str(naca2(2)),num2str(naca2(3))];
%naca3_str = [num2str(naca3(1)),num2str(naca3(2)),num2str(naca3(3))];

coo1 = fourdigit(naca1,100);
coo2 = fourdigit(naca2,100);
%coo3 = fourdigit(naca3,100);

figure(1),clf
plot(coo1(:,1),coo1(:,2),'k'),hold on
plot(coo2(:,1),coo2(:,2),'r--')
%plot(coo3(:,1),coo3(:,2),'b')
grid on,axis equal
%legend(['Caso 1 - NACA ',naca1_str],['Caso 2 - NACA ',naca2_str],['Caso 3 - NACA ',naca3_str])
%legend(['Raiz - NACA ',naca1_str],['Meio - NACA ',naca2_str],['Ponta - NACA ',naca3_str])
legend(['Caso 1 - NACA ',naca1_str],['Caso 2 - NACA ',naca2_str])
%legend(['Caso 1 - NACA ',naca1_str])
set(gca,'xlim',[0,1])
set(gca,'ylim',[-0.15,0.25])

% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)


saveas(gcf,'aoba.png')