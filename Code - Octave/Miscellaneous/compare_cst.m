clc,clear

% perfil 1 (original)
v_ex1 = [0.0100, 0.2000, 0.2000, 0.0000, 10.0000, 0.0000];
v_in1 = [0.0050, 0.0000, 0.0750, 0.0000, 5.0000, 0.0000];

% perfil 2 (otimizado)
v_ex2 = [0.0500, 0.0000, 0.2500, 0.0500, 20.0000, 0.0000];
v_in2 = [0.0300, 0.0500, 0.1000, 0.1000, 10.0000, 0.0000];

% perfil 3
v_ex3 = [0.0200, 0.2750, 0.2500, 0.0500, 13.0000, 0.0000];
v_in3 = [0.0300, 0.1000, 0.0250, -0.0000, 9.0000, 0.0000];




dat.chord = 1;
%n = dat.BPn;
%dat.n = length(v_ex1) - 2;
dat.np = 100;
dat.N1 = 0.5;
dat.N2 = 1;
dat.p_op = 1;

dat1 = dat;
dat2 = dat;
dat3 = dat;
dat1.n = length(v_ex1) - 2;
dat2.n = length(v_ex2) - 2;
dat3.n = length(v_ex3) - 2;


coo1 = run_cst_TCC2(v_ex1,v_in1,dat1);
coo2 = run_cst_TCC2(v_ex2,v_in2,dat2);
coo3 = run_cst_TCC2(v_ex3,v_in3,dat3);




figure(1),clf
plot(coo1(:,1),coo1(:,2),'k'),hold on
plot(coo2(:,1),coo2(:,2),'r--')
%plot(coo3(:,1),coo3(:,2),'b')
grid on,axis equal
%legend('Raiz','Meio','Ponta')
%legend('Caso 1','Caso 2','Caso 3')
%legend('Original','Otimizado (minimizar momento)')
%legend('Sem requisito de CM constante','Com requisito de CM constante')
legend('Caso 1','Caso 2')

set(gca,'xlim',[0,1])
set(gca,'ylim',[-0.08,0.22])

% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

saveas(gcf,'aoba1.png')


%figure(2),clf
%plot(coo1(:,1),coo1(:,2),'k--'),hold on
%plot(coo3(:,1),coo3(:,2),'r')
%grid on,axis equal
%legend('Original','Otimizado (minimizar momento)')
%
%set(gca,'xlim',[0,1])
%set(gca,'ylim',[-0.1,0.2])
%
%% Trocar separador decimal
%xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
%new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
%set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)
%
%saveas(gcf,'aoba2.png')