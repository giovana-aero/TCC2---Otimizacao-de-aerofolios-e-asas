function v = make_vector(pop,op,select2)
% Esta função transforma informações dos vetores aero e transforma em
% um vetor

v = zeros(1,length(pop));
for i = [select2]
    v(i) = pop(i).aero(op);
end

%% Remover valores referentes aos indivíduos defeituosos
%for i = [flip(select)]
%    v(i) = [];
%end
    
end