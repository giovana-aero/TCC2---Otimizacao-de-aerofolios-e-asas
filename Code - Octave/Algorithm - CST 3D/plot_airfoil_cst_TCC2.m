function plot_airfoil_cst_TCC2(coordenadas,op,loop,pop)
    
% - As duas primeiras entradas (coordenadas,op) são obrigatórias para todos os casos
% - A terceira entrada é obrigatória para op = 1 e op = 2
% - A quarta entrada é obrigatória para op = 2

if op == 0 % Ler coordenadas a partir de um arquivo de texto
    coordenadas = dlmread('coordenadas.dat');
    plot(coordenadas(:,1),coordenadas(:,2))

elseif op == 1 % Por no título apenas o número da iteração
    plot(coordenadas(:,1),coordenadas(:,2))
    title(['Iteração ' num2str(loop)])
    
elseif op == 2 % Pôr no título as informações aerodinâmicas
    plot(coordenadas(:,1),coordenadas(:,2))
    aero = pop.aero;
    title(['Iteração ' num2str(loop) ': CL = ' num2str(aero(1))...
          ', CD = ' num2str(aero(2)) ', L/D = ', num2str(aero(3)) ', CM = ' num2str(aero(4))])
    
elseif op == 3 % Fazer com estilo tracejado pra fins de comparação
    plot(coordenadas(:,1),coordenadas(:,2),'--')
   
  
elseif op == 4 % Fazer um gráfico sem nada no título
  plot(coordenadas(:,1),coordenadas(:,2))
 
end
%axis equal,grid on


% Trocar separador decimal
xl = get(gca,'XTickLabel'); yl = get(gca,'YTickLabel');
new_xl = strrep(xl(:),'.',','); new_yl = strrep(yl(:),'.',',');
set(gca,'XTickLabel',new_xl), set(gca,'YTickLabel',new_yl)

end