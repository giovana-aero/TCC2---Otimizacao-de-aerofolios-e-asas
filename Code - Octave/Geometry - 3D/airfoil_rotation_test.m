% Objetivo: ler as coordenadas de um aerofólio e rotacioná-lo (torção geométrica)
clc,clear

nome = 'coordenadas.dat'; % O formato das coordenadas deve ser o do XFOIL
th = -5; % ângulo de rotação [graus]
c = 1/4; % Centro de rotação (ao longo da linha da corda)

coo = dlmread(nome);

% Traçar o perfil original
figure(1),clf
plot(coo(:,1),coo(:,2),'k'),hold on,grid on,axis equal
c_L = [0,0;1,0];
plot(c_L(:,1),c_L(:,2),'k--') % Linha de corda

% Fazer a rotação
th = -th; % As convenções de sinal (matemática e aeronáutica) são inversas, portanto, realiza-se a troca aqui
R = [cosd(th),-sind(th);
     sind(th),cosd(th)];
coo(:,1) = coo(:,1) - c; c_L(:,1) = c_L(:,1) - c; % Transladar as coordenadas pra que a rotação ocorra em torno do ponto desejado
coo_R = (R*coo')';                  c_L_R = (R*c_L')';
coo_R(:,1) = coo_R(:,1) + c;        c_L_R(:,1) = c_L_R(:,1) + c; % Mover a geometria de volta à posição original

% Traçar o perfil rotacionado
plot(coo_R(:,1),coo_R(:,2),'b')
plot(c_L_R(:,1),c_L_R(:,2),'b--') 
scatter(c,0,'filled')

    