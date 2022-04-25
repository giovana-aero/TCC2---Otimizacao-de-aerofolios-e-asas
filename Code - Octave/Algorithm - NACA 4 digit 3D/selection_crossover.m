function n = selection_crossover(p)

r = rand*sum(p);
c = cumsum(p);
n = find(r<=c,1,'first');

end