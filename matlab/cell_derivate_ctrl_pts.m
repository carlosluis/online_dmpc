function T_ctrl_pts = cell_derivate_ctrl_pts(d)
der_mats{d} = [];
for c = d:-1:1
    aux = zeros(c, c + 1);
    for l = 1:c
       aux(l, l:l + 1) = [-c, c]; 
    end
    der_mats{c} = aux;
end

n = d; % n represents the n-th derivative of position
for k = 1:d
    aux = 1;
    for k_aux = d:-1:k
       aux = der_mats{k_aux}*aux;   
    end
    T_ctrl_pts{n}(:,:) = aux; 
    n = n - 1;
end