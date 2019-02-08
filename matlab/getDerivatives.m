function der_matrices = getDerivatives(deg)
der_matrices{deg} = [];
for c = deg:-1:1
    aux = zeros(c,c+1);
    for l = 1:c
       aux(l,l:l+1) = [-c c]; 
    end
    der_matrices{c} = aux;
end