
% Nota: o vetor select não é utilizado aqui. resolver isso depois




function pop = fitness_naca4_3D(pop,dat,select2)
    
for P = 1:dat.cases % Fazer loops para cada condição de voo
    % CL (aero(1))
    if dat.coeff_op(P,1) ~= '!'
        F1 = dat.coeff_F(P,1);
        if dat.coeff_op(P,1) == '^' % Maior CL possível
            temp = make_vector_TCC2_3D(pop,P,1,1,select2); % Vetor com todos os CLs
            for i = [select2]
                pop(i).score = pop(i).aero(P,1)/max(temp)*F1 + pop(i).score;
            end
        elseif dat.coeff_op(P,1) == 'o' % Alcançar um CL específico
            CLtgt = dat.coeff_val(P,1);
%            if CLtgt == 0 && dat.aoa(P) == 0 % Se o CL alvo for igual a zero e o ângulo de ataque for nulo 
%                for i = [select2]
%                    pop(i).score = 1 + pop(i).score;
%                end
            if CLtgt == 0 % Se o CL alvo for igual a zero
                temp = make_vector_TCC2_3D(pop,P,1,1,select2); % Vetor com todos os CLs
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
            
            
            % ATENÇÃO
            
            
        elseif dat.coeff_op(P,1) == 'q' % Alcançar uma força de sustentação específica
            Ltgt = dat.coeff_val(P,1);
            if Ltgt == 0 % Se o L alvo for igual a zero   *1/2*dat.rho(P)*dat.v_ref(P)^2*pop.S
                temp = make_vector_TCC2_3D(pop,P,1,2,select2,dat.v_ref,dat.rho); % Vetor com todos os Ls
                for i = [select2]
                    pop(i).score = (1-abs(temp(i))/max(abs(temp)))*F1 + pop(i).score;
                end 
            elseif Ltgt > 0 % Se o L alvo for positivo
                for i = [select2]
                    L_pop = pop(i).aero(P,1)*1/2*dat.rho(P)*dat.v_ref(P)^2*pop(i).S;
                    if pop(i).aero(P,1) >= Ltgt % Se o L for maior/igual que o L alvo
                        pop(i).score = Ltgt/L_pop*F1 + pop(i).score;
                    else % Se o L for menor que o L alvo
                        pop(i).score = L_pop/Ltgt*F1 + pop(i).score;
                    end
                end
            else % Se o L alvo for negativo
                for i = [select2]
                    L_pop = pop(i).aero(P,1)*1/2*dat.rho(P)*dat.v_ref(P)^2*pop(i).S;
                    if pop(i).aero(P,1) <= Ltgt % Se o L for menor/igual que o L alvo
                        pop(i).score = Ltgt/L_pop*F1 + pop(i).score;
                    else % Se o L for maior que o L alvo
                        pop(i).score = L_pop/Ltgt*F1 + pop(i).score;
                    end
                end
            end
            
            
            
        elseif dat.coeff_op(P,1) == '#' % Maior sustentação possível
            temp = make_vector_TCC2_3D(pop,P,1,2,select2,dat.v_ref,dat.rho); % Vetor com todos os Ls
            for i = [select2]
                pop(i).score = temp(i)/max(temp)*F1 + pop(i).score;
            end

    

            % ATENÇÃO
        
    
    
            
        end
    end
    % CD (aero(2))
    if dat.coeff_op(P,2) ~= '!'
        F2 = dat.coeff_F(P,2);
        if dat.coeff_op(P,2) == '^' % Menor CD possível (tende a zero, mas nunca igual a ou menor que zero)
            temp = make_vector_TCC2_3D(pop,P,2,1,select2);
            for i = [select2]
                if pop(i).aero(P,2) <= 0
                    pop(i).score = 0; % Punir asas com resultados irreais
                else
                    pop(i).score = (1-pop(i).aero(P,2)/max(temp))*F2 + pop(i).score;
                end
            end
        elseif dat.coeff_op(P,2) == 'o' % Alcançar um CD específico (apenas valores positivos)
            CDtgt = dat.coeff_val(P,2);
            for i = [select2]
                if pop(i).aero(P,2) >= CDtgt % Se o CD for maior/igual que o CD alvo
                    pop(i).score = CDtgt/pop(i).aero(P,2)*F2 + pop(i).score;
                else % Se o CD for menor que o CD alvo
                    pop(i).score = pop(i).aero(P,2)/CDtgt*F2 + pop(i).score;
                end
            end
            
            
            % ATENÇÃO
            
            
        elseif dat.coeff_op(P,2) == 'q' % Alcançar uma força de arrasto específica 
            Dtgt = dat.coeff_val(P,2);
            for i = [select2]
                D_pop = pop(i).aero(P,2)*1/2*dat.rho(P)*dat.v_ref(P)^2*pop(i).S;
                if D_pop >= Dtgt % Se o D for maior/igual que o D alvo
                    pop(i).score = Dtgt/D_pop*F1 + pop(i).score;
                else % Se o D for menor que o D alvo
                    pop(i).score = D_pop/Dtgt*F1 + pop(i).score;
                end
            end
        
            
        elseif dat.coeff_op(P,2) == '#' % Menor arrasto possível
            temp = make_vector_TCC2_3D(pop,P,2,2,select2,dat.v_ref,dat.rho);
            for i = [select2]
                if pop(i).aero(P,2) <= 0
                    pop(i).score = 0; % Punir asas com resultados irreais
                else
                    pop(i).score = (1-temp(i)/max(temp))*F2 + pop(i).score;
                end
            end
            
            
            
            
            % ATENÇÃO
            
            
        end
    end
    % L/D (aero(3))
    if dat.coeff_op(P,3) ~= '!'
        F3 = dat.coeff_F(P,3);
        if dat.coeff_op(P,3) == '^' % Maior L/D possível
            temp = make_vector_TCC2_3D(pop,P,3,1,select2); % Vetor com todos os L/Ds
            for i = [select2]
                pop(i).score = pop(i).aero(P,3)/max(temp)*F3 + pop(i).score;
            end
        elseif dat.coeff_op(P,3) == 'o' % Alcançar um L/D específico
            LDtgt = dat.coeff_val(P,3);
%            if LDtgt == 0 && dat.aoa(P) == 0 % Se o LD alvo for igual a zero e o ângulo de ataque for nulo 
%                for i = [select2]
%                    pop(i).score = 1 + pop(i).score;
%                end
            if LDtgt == 0 % Se o L/D alvo for igual a zero
                temp = make_vector_TCC2_3D(pop,P,3,1,select2); % Vetor com todos os L/Ds
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
        if dat.coeff_op(P,4) == 'o' % Alcançar um CM específico
            CMtgt = dat.coeff_val(P,4);
    %        if CMtgt == 0 && dat.aoa(P) == 0 % Se o CM alvo for igual a zero e o ângulo de ataque for nulo
    %            for i = [select2]
    %                pop(i).score = 1 + pop(i).score;
    %            end
            if CMtgt == 0 % Se o CM alvo for igual a zero
                temp = make_vector_TCC2_3D(pop,P,4,1,select2); % Vetor com todos os CMs
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
    
    
        % ATENÇÃO
        
    
        elseif dat.coeff_op(P,4) == 'q' % Alcançar um M específico
            Mtgt = dat.coeff_val(P,4);
    %        if CMtgt == 0 && dat.aoa(P) == 0 % Se o CM alvo for igual a zero e o ângulo de ataque for nulo
    %            for i = [select2]
    %                pop(i).score = 1 + pop(i).score;
    %            end
            if Mtgt == 0 % Se o M alvo for igual a zero
                temp = make_vector_TCC2_3D(pop,P,4,3,select2,dat.v_ref,dat.rho); % Vetor com todos os Ms
                for i = [select2]
                    pop(i).score = (1-abs(temp(i))/max(abs(temp)))*F4 + pop(i).score;
                end
            elseif Mtgt > 0 % Se o M alvo for positivo
                for i = [select2]
                    M_pop = pop(i).aero(P,4)*1/2*dat.rho(P)*dat.v_ref(P)^2*pop(i).S*pop(i).mac;
                    if pop(i).aero(P,4) >= Mtgt % Se o M for maior/igual que o M alvo
                        pop(i).score = Mtgt/M_pop*F4 + pop(i).score;
                    else % Se o M for menor que o M alvo
                        pop(i).score = M_pop/Mtgt*F4 + pop(i).score;
                    end
                end
            else % Se o M alvo for negativo
                for i = [select2]
                    M_pop = pop(i).aero(P,4)*1/2*dat.rho(P)*dat.v_ref(P)^2*pop(i).S*pop(i).mac;
                    if pop(i).aero(P,4) <= Mtgt % Se o M for menor/igual que o M alvo
                        pop(i).score = Mtgt/M_pop*F4 + pop(i).score;
                    else % Se o M for maior que o M alvo
                        pop(i).score = M_pop/Mtgt*F4 + pop(i).score;
                    end
                end
            end
        
            
            
            % ATENÇÃO
        
        
        end
    end
end

% Funções objetivas c e k para o coeficiente de momento
if dat.coeff_op(1,4) == 'c' || dat.coeff_op(1,4) == 'k' 
    for i = [select2]
        temp = 0;
        if dat.coeff_op(1,4) == 'c' % Se o valor for arbitrário
            CMtgt = mean(pop(i).aero(:,4));
        else dat.coeff_op(1,4) == 'k' % Se o valor for específico
            CMtgt = dat.coeff_val(1,4);
        end
        for P = 1:dat.cases
            if CMtgt == 0 % Se o CM alvo for igual a zero
                temp2 = make_vector_TCC2_3D(pop,P,4,1,select2); % Vetor com todos os CMs
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