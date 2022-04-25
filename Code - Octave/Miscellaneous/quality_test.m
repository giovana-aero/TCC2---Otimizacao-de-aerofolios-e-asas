clc,clear

v_ex1 = [0.01,0.1,0.1,0.1,10,0];
v_in1 = [0.01,0.2,0.2,0.1,10,0];

% pegar informações do struct
dat.chord = 1;
%n = dat.BPn;
n = length(v_ex1) - 2;
dat.np = 100;
dat.N1 = 0.5;
dat.N2 = 1;
dat.p_op = 1;

[coo1,mean] = run_cst_TCC2(v_ex1,v_in1,dat,0);
check = quality(coo1,dat);disp(check)

figure(1),clf
plot(coo1(:,1),coo1(:,2),'k'),axis equal,grid on,hold on
plot(mean(:,1),mean(:,2),'b')

set(gca,'xlim',[0,1])
set(gca,'ylim',[-0.2,0.2])

% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

saveas(gcf,'aoba.png')

