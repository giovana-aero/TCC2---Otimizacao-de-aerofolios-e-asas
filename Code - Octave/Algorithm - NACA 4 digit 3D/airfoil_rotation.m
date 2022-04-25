function coo_R = airfoil_rotation(coo,th)

% Centro de rotação (sempre no ponto de quarto de corda como padrão)
c = 1/4;

% Fazer a rotação
th = -th; % As convenções de sinal (matemática e aeronáutica) são inversas, portanto, realiza-se a troca aqui
R = [cosd(th),-sind(th);
     sind(th),cosd(th)];
coo(:,1) = coo(:,1) - c;  % Transladar as coordenadas pra que a rotação ocorra em torno do ponto desejado
coo_R = (R*coo')';                  
coo_R(:,1) = coo_R(:,1) + c;     % Mover a geometria de volta à posição original

end