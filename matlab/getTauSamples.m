function Tau = getTauSamples(T_segment,t_vec,dim, l)
% Tau will return a 3D matrix, where the last subscript represents from
% which segment of the Bezier curve the sample is to be taken
T_final = l*T_segment;
Tau{l} = [];
for k = 1:length(t_vec)
    if t_vec(k) == T_final
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
end