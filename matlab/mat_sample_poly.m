function T_sample_poly = mat_sample_poly(T_segment, t_vec, dim, l)
% Create a cell array Tau with the appropriate samples 
% to be taken from each segment
T_final = l*T_segment;
Tau{l} = [];
for k = 1:length(t_vec)
    if abs(t_vec(k) - T_final) < 1e-4
       segment = l;
       t = T_segment;
    else
        segment = floor(t_vec(k)/T_segment) + 1;
        t = rem(t_vec(k),T_segment);
    end
    row = size(Tau{segment},1) + 1;
    for n = 0:dim
        Tau{segment}(row,n+1) = (t/T_segment)^n;
    end
end

T_sample_poly = get_sample_mat(Tau,length(t_vec),dim);
end