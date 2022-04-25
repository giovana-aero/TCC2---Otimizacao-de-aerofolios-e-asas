function n = selection_crossover(p)

r = rand*sum(p);
c = cumsum(p);
n = find(r<=c,1,'first');

end


%weights = [0.12,0.04,0.16,0.08,0.36,0.24];
%c =    0.1200   0.1600   0.3200   0.4000   0.7600   1.0000;