function coo_intp = airfoil_interpolation(coo_r,coo_t,op,d,tw)
% Reescrita do script test_interpolation_v02 
% Adaptação para ser utilizada no script de malha do APAME:
% - Todos os aerofólios têm coordenadas na mesma configuração/tamanho, portanto, 
%   não é necessário utilizar as splines

% Nota: o tamanho das coordenadas é igual em dois sentidos:
% - Número total de nós
% - Número de nós para cada superfície (extradorso e intradorso). Isso faz com 
%   que o bordo de ataque sempre fique na mesma posição


% Encontrar o bordo de ataque de ambos os aerofólios
delta_r = zeros(1,size(coo_r,1)-1);
for i = 1:(size(coo_r,1)-1)
    delta_r(i) = coo_r(i+1,1) - coo_r(i,1);
end
delta_b_r = delta_r >= 0;
LE = find(delta_b_r == 1,1);

% Separar superfícies (ordenadas)
ex_r = coo_r(1:LE,:);
in_r = coo_r(LE:end,:);
ex_t = coo_t(1:LE,:);
in_t = coo_t(LE:end,:);

% Fazer interpolações
% 0 -> raiz
% 1 -> ponta
ex_intp = (ex_r*(1-op) + ex_t*op);
in_intp = (in_r*(1-op) + in_t*op);

% Montar coordenadas
coo_intp = [ex_intp;in_intp(2:end,:)];

% Fazer rotação do perfil
coo_intp = airfoil_rotation(coo_intp,tw);

if nargin > 3 % Passar as ordenadas pro eixo z e aplicar um deslocamento ao longo do eixo y
%    n = size(coo_intp,1);
    coo_intp = [coo_intp(:,1),repmat(d,size(coo_intp,1),1),coo_intp(:,2)];
    
%NODE(1:sec_N,:) = [coo_t(:,1),repmat(-b/2,sec_N,1),coo_t(:,2)];
else % Retornar as coordenadas no formato 2D comum
    

end




