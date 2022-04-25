function dat = error_check_cst_3D(dat)

% Forçar a variável chord para ter valor unitário. Isso é importante para evitar
% conflito com as especificações de tamanhos de corda da planta da asa
dat.chord = 1;

% Parâmetros do algoritmo: tamanho da população
if rem(dat.N,2) ~= 0
    warning('Tamanho da população deve ser par. O valor inserido será alterado'),pause(5)
    dat.N = dat.N + 1;
end

% Parâmetros do algoritmo: chance de mutação
if dat.mu < 0 || dat.mu > 1
    error('Chance de mutação deve estar no intervalo 0 <= dat.mu <= 1')
end

% Parâmetros do algoritmo: número de iterações
if dat.iter < 1 || dat.iter ~= floor(dat.iter)
    error('Número de iterações do algoritmo deve ser positivo, não nulo e inteiro')
end

% Parâmetros da geometria: proporção de tipos de planta (exclusivo ao algoritmo geral)
if isfield(dat,'planf_op')
    if dat.planf_op < 0 || dat.planf_op > 1
        error('A proporção de tipos de planta deve estar no intervalo 0 <= dat.planf_op <= 1')
    end
end

% Parâmetros da geometria: configuração 'L' da corda do meio e do perfil do meio
if dat.planf_op ~= 0 && ismember('L',dat.c_m_ext_in) && ismember('L',dat.m_ext_m)
	error("Configurações da geometria transformarão todas as asas bitrapezoidais em trapezoidais simples (rever opções 'L')")
end

% Parâmetros da geometria: extensão da envergadura da primeira seção b1
if dat.b1_ext_in(1) >= dat.b_ext_in(2) && dat.planf_op ~= 0
	error('O valor mínimo da extensão de b1 deve sempre ser menor que o valor máximo da extensão de b')
end

% Parâmetros da geometria: extensão da envergadura da primeira seção b1
if dat.b1_ext_in(3) <= 0 && dat.planf_op ~= 0
	error('A separação mínima entre b1 e b deve ser um valor positivo não nulo')
end

% Parâmetros da geometria: extensão da corda do meio
if dat.c_m_ext_in(2) > dat.c_r_ext_in(2) && dat.planf_op ~= 0
	error('O valor máximo da extensão de c_m deve ser menor ou igual ao valor da extensão máxima da extensão de c_r')
end

% Parâmetros da geometria: extensão da corda do meio
if dat.c_m_ext_in(1) < dat.c_t_ext_in(1) && dat.planf_op ~= 0
	error('O valor mínimo da extensão de c_m deve ser maior ou igual ao valor mínimo da extensão de c_t')
end

% Erros relacionados aos aerofólios?
if dat.symm_op_r < 0 || dat.symm_op_r > 1 || dat.symm_op_m < 0 || dat.symm_op_m > 1 || dat.symm_op_t < 0 || dat.symm_op_t > 1
	error('Opção de perfis simétricos deve estar definida no intervalo 0 <= dat.symm_op <= 1')
end

% Parâmetros das simulações: número de condições de voo
if dat.cases < 1 || dat.cases ~= floor(dat.cases)
	error('Número de condições de voo deve ser positivo, não nulo e inteiro')
end

% Parâmetros das simulações: velocidades de referência
if length(dat.v_ref) < dat.cases
	error('Não há velocidades de referência suficientes para todas as condições de voo')
end

% Parâmetros das simulações: densidades do ar
if length(dat.rho) < dat.cases
	error('Não há densidades do ar suficientes para todas as condições de voo')
end 

% Parâmetros das simulações: pressões atmosféricas
if length(dat.p_atm) < dat.cases
	error('Não há pressões atmosféricas suficientes para todas as condições de voo')
end

% Parâmetros das simulações: números de Mach
if length(dat.mach) < dat.cases
	error('Não há números de Mach suficientes para todas as condições de voo')
end

% Parâmetros das simulações: números de Reynolds
if length(dat.reynolds) < dat.cases 
    error('Não há números de Reynolds suficientes para todas as condições de voo')
end

% Parâmetros das simulações: ângulos de ataque
if length(dat.aoa) < dat.cases 
    error('Não há ângulos de ataque suficientes para todas as condições de voo')
end

%if length(dat.iter_sim) < dat.cases
%    error('Não há números de iterações (XFOIL) suficientes para todas as condições de voo')
%end

% % Trocar valores nulos e negativos de números de iteração pelo valor padrão do XFOIL
% for i = 1:dat.cases
	% if dat.iter_sim(i) == 0
		% dat.iter_sim(i) = 10;
	% end
% end

% Parâmetros das simulações: configuração das funções objetivas
if size(dat.coeff_op,1) < dat.cases || size(dat.coeff_val,1) < dat.cases || size(dat.coeff_F,1) < dat.cases
    error('Configurações das funções objetivas são insuficientes para todas as condições de voo')
end

% Parâmetros das simulações: função objetiva 'q' para força de sustentação, força de arrasto e momento de arfagem
if ismember('q',dat.coeff_op(:,3))
	error('Função objetiva q vale apenas para a sustentação, arrasto e momento')
end

% parâmetros das simulações: função objetiva '|' para força de sustentação e de arrasto
if ismember('#',dat.coeff_op(:,3:4))
	error('Função objetiva # vale apenas para sustentação e arrasto')
end

% Parâmetros das simulações: funções objetivas de coeficiente de momento
if ismember('c',dat.coeff_op(:,1:3)) || ismember('k',dat.coeff_op(:,1:3))
    error('Funções objetivas c e k valem apenas para o coeficiente de momento')
end

% Parâmetros das simulações: funções objetivas de coeficientes de momento referentes a múltiplas condições de voo
if dat.cases == 1 && dat.coeff_op(1,4) == 'c' || dat.cases == 1 && dat.coeff_op(1,4) == 'k'
    warning('Função objetiva c/k vale apenas para múltiplas condições de voo. A mesma será ignorada.')
    dat.coeff_op(1,4) = '!'; pause(5)
end

if sum(sum(dat.coeff_op(1:dat.cases,:) == '!')) == numel(dat.coeff_op(1:dat.cases,:)),
    error('Ao menos uma das funções objetivas deve estar ativa')
end

% PArâmetros das simulações: funções objetivas (geral)
T = dat.coeff_op;
for i = 1:dat.cases
	for j = 1:4
		if T(i,j) ~= '!' && T(i,j) ~= '^' && T(i,j) ~= 'o' && T(i,j) ~= 'q' && T(i,j) ~= '#' && T(i,j) ~= 'c' && T(i,j) ~= 'k'
			error('A matriz dat.coeff_op deve ser definida com opções !, ^, o, q, #, c ou k')
		end
	end
end

if dat.cases > 1
	for i = 2:dat.cases
        if sum(dat.coeff_op(i,:) == '!') == 4 && dat.coeff_op(1,4) ~= 'c' && dat.coeff_op(1,4) ~= 'k'
			error(['Condição de voo ' num2str(i) ' não tem função objetiva definida'])
		end
	end
end

for P = 1:dat.cases
    if dat.coeff_op(P,2) == 'o' && dat.coeff_val(P,2) < 0
        error(['Condição de voo ' num2str(P) ': CD alvo deve ser maior que zero'])
    elseif dat.coeff_op(P,2) == 'o' && dat.coeff_val(P,2) == 0
        dat.coeff_op(P,2) = '^';
        warning(['Condição de voo ' num2str(P) ': CD = 0 - Função objetiva de arrasto trocada de o para ^']),pause(5)
    end
end

if dat.cases > 1 && dat.coeff_op(1,4) == 'c' || dat.cases > 1 && dat.coeff_op(1,4) == 'k'
	if sum(dat.coeff_op(2:end,4) ~= '!') ~= 0
		warning('Função objetiva de momento c/k: funções objetivas das condições de voo subsequentes serão ignoradas'),pause(5)
		dat.coeff_op(2:dat.cases,4) = repmat('!',dat.cases-1,1);
	end	
end
