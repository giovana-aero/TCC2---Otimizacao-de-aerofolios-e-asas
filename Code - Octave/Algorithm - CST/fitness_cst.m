function pop = fitness_cst(pop,dat,select2)
    
for P = 1:dat.cases % Fazer loops para cada condição de voo
    % CL (aero(1))
    if dat.coeff_op(P,1) ~= '!'
        F1 = dat.coeff_F(P,1);
        if dat.coeff_op(P,1) == '^' % Maior CL possível'
            temp = make_vector_TCC2(pop,P,1,select2); % Vetor com todos os CLs
            for i = [select2]
                pop(i).score = pop(i).aero(P,1)/max(temp)*F1 + pop(i).score;
            end
        elseif dat.coeff_op(P,1) % Alcançar um CL específico
            CLtgt = dat.coeff_val(P,1);
            if CLtgt == 0 && dat.aoa(P) == 0 && dat.symm_op == 1 % Se o CL alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
                for i = [select2]
                    pop(i).score = 1 + pop(i).score;
                end
            elseif CLtgt == 0 % Se o CL alvo for igual a zero
                temp = make_vector_TCC2(pop,P,1,select2); % Vetor com todos os CLs
                for i = [select2]
                    pop(i).score = (1-abs(pop(i).aero(P,1))/max(abs(temp)))*F1  + pop(i).score;
                end 
            elseif CLtgt > 0 % Se o CL alvo for positivo
                for i = [select2]
                    if pop(i).aero(P,1) >= CLtgt % Se o CL for maior/igual que o CL alvo
                        pop(i).score = CLtgt/pop(i).aero(P,1)*F1 + pop(i).score;
                    else % Se o CL for menor que o CL alvo
                        pop(i).score = pop(i).aero(P,1)/CLtgt*F1 + pop(i).score;
                    end
                end
            else % Se o CL alvo for negativo
                for i = [select2]
                    if pop(i).aero(P,1) <= CLtgt % Se o CL for menor/igual que o CL alvo
                        pop(i).score = CLtgt/pop(i).aero(P,1)*F1 + pop(i).score;
                    else % Se o CL for maior que o CL alvo
                        pop(i).score = pop(i).aero(P,1)/CLtgt*F1 + pop(i).score;
                    end
                end
            end
        end
    end
    % CD (aero(2))
    if dat.coeff_op(P,2) ~= '!'
        F2 = dat.coeff_F(P,2);
        if dat.coeff_op(P,2) == '^' % Menor CD possível (tende a zero, mas nunca igual a ou menor que zero)
            temp = make_vector_TCC2(pop,P,2,select2);
            for i = [select2]
                if pop(i).aero(P,2) <= 0
                    pop(i).score = 0; % Punir aerofólios com resultados irreais
                else
                    pop(i).score = (1-pop(i).aero(P,2)/max(temp))*F2 + pop(i).score;
                end
            end
        elseif dat.coeff_op(P,2) % Alcançar um CD específico (apenas valores positivos)
            CDtgt = dat.coeff_val(P,2);
            for i = [select2]
                if pop(i).aero(P,2) >= CDtgt % Se o CD for maior/igual que o CD alvo
                    pop(i).score = CDtgt/pop(i).aero(P,2)*F2 + pop(i).score;
                else % Se o CD for menor que o CD alvo
                    pop(i).score = pop(i).aero(P,2)/CDtgt*F2 + pop(i).score;
                end
            end
        end
    end
    % L/D (aero(3))
    if dat.coeff_op(P,3) ~= '!'
        F3 = dat.coeff_F(P,3);
        if dat.coeff_op(P,3) == '^' % Maior L/D possível
            temp = make_vector_TCC2(pop,P,3,select2); % Vetor com todos os L/Ds
            for i = [select2]
                pop(i).score = pop(i).aero(P,3)/max(temp)*F3 + pop(i).score;
            end
        elseif dat.coeff_op(P,3) == 'o' % Alcançar um L/D específico
            LDtgt = dat.coeff_val(P,3);
            if LDtgt == 0 && dat.aoa(P) == 0 && dat.symm_op == 1 % Se o LD alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
                for i = [select2]
                    pop(i).score = 1 + pop(i).score;
                end
            elseif LDtgt == 0 % Se o L/D alvo for igual a zero
                temp = make_vector_TCC2(pop,P,3,select2); % Vetor com todos os L/Ds
                for i = [select2]
                    pop(i).score = (1-abs(pop(i).aero(P,3))/max(abs(temp)))*F3 + pop(i).score;
                end
            elseif LDtgt > 0 % Se o L/D alvo for positivo
                for i = [select2]
                    if pop(i).aero(P,3) >= LDtgt % Se o L/D for maior/igual que o L/D alvo
                        pop(i).score = LDtgt/pop(i).aero(P,3)*F3 + pop(i).score;
                    else % Se o L/D for menor que o L/D alvo
                        pop(i).score = pop(i).aero(P,3)/LDtgt*F3 + pop(i).score;
                    end
                end
            else % Se o L/D alvo for negativo
                if pop(i).aero(P,3) <= LDtgt % Se o L/D for menor/igual que o L/D alvo
                    pop(i).score = LDtgt/pop(i).aero(P,3)*F3 + pop(i).score;
                else % Se o L/D for maior que o L/D alvo
                    pop(i).score = pop(i).aero(P,3)/LDtgt*F3 + pop(i).score;
                end
            end
        end
    end
    % CM (aero(4))
    if dat.coeff_op(P,4) ~= '!' && dat.coeff_op(P,4) ~= 'c' && dat.coeff_op(P,4) ~= 'k'
        F4 = dat.coeff_F(P,4);
        CMtgt = dat.coeff_val(P,4);
        if CMtgt == 0 && dat.aoa(P) == 0 && dat.symm_op == 1 % Se o CM alvo for igual a zero, o ângulo de ataque for nulo e todos os perfis forem simétricos
            for i = [select2]
                pop(i).score = 1 + pop(i).score;
            end
        elseif CMtgt == 0 % Se o CM alvo for igual a zero
            temp = make_vector_TCC2(pop,P,4,select2); % Vetor com todos os CMs
            for i = [select2]
                pop(i).score = (1-abs(pop(i).aero(P,4))/max(abs(temp)))*F4 + pop(i).score;
            end
        elseif CMtgt > 0 % Se o CM alvo for positivo
            for i = [select2]
                if pop(i).aero(P,4) >= CMtgt % Se o CM for maior/igual que o CM alvo
                    pop(i).score = CMtgt/pop(i).aero(P,4)*F4 + pop(i).score;
                else % Se o CM for menor que o CM alvo
                    pop(i).score = pop(i).aero(P,4)/CMtgt*F4 + pop(i).score;
                end
            end
        else % Se o CM alvo for negativo
            for i = [select2]	
                if pop(i).aero(P,4) <= CMtgt % Se o CM for menor/igual que o CM alvo
                    pop(i).score = CMtgt/pop(i).aero(P,4)*F4 + pop(i).score;
                else % Se o CM for maior que o CM alvo
                    pop(i).score = pop(i).aero(P,4)/CMtgt*F4 + pop(i).score;
                end
            end
        end
    end
end

% Funções objetivo c e k para o coeficiente de momento
if dat.coeff_op(1,4) == 'c' || dat.coeff_op(1,4) == 'k' 
    for i = [select2]
        temp = 0;
        if dat.coeff_op(1,4) == 'c' % Se o valor for arbitrário
            CMtgt = mean(pop(i).aero(:,4));
        else % Se o valor for específico
            CMtgt = dat.coeff_val(1,4);
        end
        for P = 1:dat.cases
            if CMtgt == 0 % Se o CM alvo for igual a zero
                temp2 = make_vector_TCC2(pop,P,4,select2); % Vetor com todos os CMs
                temp = (1-abs(pop(i).aero(P,4))/max(abs(temp2))) + temp;
            elseif CMtgt > 0 % Se o CM alvo for positivo
                if pop(i).aero(P,4) >= CMtgt % Se o CM for maior/igual que o CM alvo
                    temp = CMtgt/pop(i).aero(P,4) + temp;
                else % Se o CM for menor que o CM alvo
                    temp = pop(i).aero(P,4)/CMtgt + temp;
                end
            else % Se o CM alvo for negativo
                if pop(i).aero(P,4) <= CMtgt % Se o CM for menor/igual que o CM alvo
                    temp = CMtgt/pop(i).aero(P,4) + temp;
                else % Se o CM for maior que o CM alvo
                    temp = pop(i).aero(P,4)/CMtgt + temp;
                end
            end
        end            
        pop(i).score = (temp/dat.cases)*dat.coeff_F(1,4) + pop(i).score; 
    end
end