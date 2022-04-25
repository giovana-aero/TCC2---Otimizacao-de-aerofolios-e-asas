clc,clear

% asa 1
pop1.type = 1;
pop1.b = 12.;
pop1.b1 = 8;
pop1.c_r = 1.8;
pop1.c_m = 1.4;
pop1.c_t = 1.4;
pop1.sweep = [];
pop1.sweep1 = 4.;
pop1.sweep2 = 11;
pop1.tw_m = 'L';
pop1.tw_t = -3.5;

% asa 2
pop2.type = 0;
pop2.b = 14;
pop2.b1 = 10.5;
pop2.c_r = 1.;
pop2.c_m = 1.;
pop2.c_t = 0.5;
pop2.sweep = 8;
pop2.sweep1 = 0;
pop2.sweep2 = 8;
pop2.tw_m = 'L';
pop2.tw_t = -1.75;

% asa 3
pop3.type = 0;
pop3.b = 14;
pop3.b1 = 10.5;
pop3.c_r = 1.;
pop3.c_m = 1.;
pop3.c_t = 0.5;
pop3.sweep = 8;
pop3.sweep1 = 0;
pop3.sweep2 = 8;
pop3.tw_m = 'L';
pop3.tw_t = -1.75;



% gráficos
figure(1),clf
plot_planform(pop1,0,'k'),hold on,grid on,axis equal
plot_planform(pop2,0,'r--')
%plot_planform(pop3,0,'b')
%legend('Caso 1','Caso 2','Caso 3')
legend('Caso 1','Caso 2')



set(gca,'xlim',[-7.5,7.5])
set(gca,'ylim',[-3,2])

% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

saveas(gcf,'aoba.png')