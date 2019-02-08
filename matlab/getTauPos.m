function Tau = getTauPos(h,T_segment,dim)
t = 0:h:T_segment;
for k =1:length(t)
    for n = 0:dim
        Tau(k,n+1) = (t(k)/T_segment)^n;
    end
end
end