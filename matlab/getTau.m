function Tau = getTau(h,t_s)
t = 0:h:t_s;
for k =1:length(t)
    Tau(k,:) = [1 t(k) t(k)^2 t(k)^3 t(k)^4];
end
end