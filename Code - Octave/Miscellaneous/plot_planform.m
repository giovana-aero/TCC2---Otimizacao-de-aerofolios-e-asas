function plot_planform(pop,op,color)

% op = 1: desenhar apenas o lado direito, op ~= 1:  desenhar ambos os lados

c_r = pop.c_r;
c_t = pop.c_t;
tw_t = pop.tw_t;
b = pop.b;

if nargin == 1
    op = 1;
    color = 'k';
elseif nargin == 2
    color = 'k';
end

axis equal,grid on
if pop.type == 0 % Se a asa for trapezoidal simples
    
    if pop.sweep == 'Z'
        x = [0,(c_r - c_t*cosd(tw_t))/2,(c_r + c_t*cosd(tw_t))/2,c_r];
    else
        x = [0,b/2*tand(pop.sweep),b/2*tand(pop.sweep)+c_t*cosd(tw_t),c_r];
    end
    y = [0,b/2,b/2,0];
    if op ~= 1 % Fazer o outro lado
        x = [x,flip(x(1:end-1))];
        y = [y,-1*flip(y(1:end-1))];
    end
    plot(y,-x,color)
    
    

%    % Bordo  de ataque da raiz ao bordo de ataque da ponta
%    plot([0,(c_r - c_t*cosd(tw_t))/2],[0,b/2],color)
%    % Bordo de ataque da ponta ao bordo de fuga da ponta
%    plot([(c_r - c_t*cosd(tw_t))/2,(c_r + c_t*cosd(tw_t))/2],[b/2,b/2],color)
%    % Bordo de fuga da ponta ao bordo de fuga da raiz
%    plot([(c_r + c_t*cosd(tw_t))/2,c_r],[b/2,0],color)


elseif pop.type == 1 % Se a asa for trapezoidal dupla
    
    b1 = pop.b1;
    c_m = pop.c_m;
    tw_m = pop.tw_m;
    % Considerar opção 'L' da corda do meio
    if c_m == 'L'
        c_m = 2*(c_t-c_r)/b*b1/2 + c_r;
    end
    % Considerar opção 'L' da torção do meio
    if tw_m == 'L' && ~ischar(tw_t) % Definir a torção do meio em termos da torção da ponta
        tw_m = tw_t*b1/b;
    elseif tw_t == 'L' && ~ischar(tw_m) % Definir a torção da ponta em termos da torção do meio
        tw_t = tw_m*b/b1;
    end
    
    
    if pop.sweep1 == 'Z' % Fazer com que o enflechamento da linha c/2 seja sempre zero
        dx_m = (c_r - c_m*cosd(tw_m))/2; % Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
    else % Usar enflechamento especificado pelo usuário
        dx_m = b1/2*tand(pop.sweep1);
    end
    if pop.sweep2 == 'Z' % Fazer com que o enflechamento da linha c/2 seja sempre zero
        dx_t = (c_m - c_t*cosd(tw_t))/2 + dx_m; % Ligeiro translado pra quando as cordas forem diferentes (é preciso considerar a torção na ponta)
%        sweep2_angle = atand((c_t-c_m)/(b2/2));
    else % Usar enflechamento especificado pelo usuário
        b2 = b - b1;
        dx_t = b2/2*tand(pop.sweep2) + dx_m;
%        sweep2_angle = sweep2;
    end
    
    x = [0,dx_m,dx_t,dx_t+c_t*cosd(tw_t),dx_m+c_m*cosd(tw_m),c_r];
%    if pop.sweep1 == 'Z' && pop.sweep2 == 'Z'
%        x = [0,(c_r - c_m*cosd(tw_m))/2,(c_r - c_t*cosd(tw_t))/2,(c_r + c_t*cosd(tw_t))/2,(c_r + c_m*cosd(tw_m))/2,c_r];
%    elseif pop.sweep1 == 'Z' && pop.sweep2 ~= 'Z'
%        x = [0,(c_r - c_m*cosd(tw_m))/2,b2/2*tand(sweep2)+dx_m,b2/2*tand(sweep2)+dx_m+c_t*cosd(tw_t),(c_r + c_m*cosd(tw_m))/2,c_r];
%    elseif pop.sweep1 ~= 'Z' && pop.sweep2 == 'Z'
%        
%    elseif pop.sweep1 == 'Z' && pop.sweep2 ~= 'Z'
%        
%    end
    y = [0,b1/2,b/2,b/2,b1/2,0];
    if op ~= 1 % Fazer o outro lado
        x = [x,flip(x(1:end-1))];
        y = [y,-1*flip(y(1:end-1))];
    end
    plot(y,-x,color)


end