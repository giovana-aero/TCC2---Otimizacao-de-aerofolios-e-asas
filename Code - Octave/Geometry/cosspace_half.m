function x = cosspace_half(startP,endP,N)
    
% Versão modificada da função cosspace. Gera apenas metade da distribuição.


x = zeros(1,N); x(N) = endP;
angleInc = pi/(N-1)/2;
    
curAngle=angleInc;
for i = 2:N-1
    x(i)=endP*(1-cos(curAngle));
    curAngle=curAngle+angleInc;
end
     
end