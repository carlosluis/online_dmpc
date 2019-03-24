function pass = reached_goal(p,pf,error_tol,N)

if(N>1)
    differ = squeeze(p(:,end,:))- squeeze(pf);
    max_dist = max(sqrt(sum(differ.^2,1)));
else
    differ = p(:,end) - pf';
    max_dist = max(sqrt(sum(differ.^2,1)));
end
    
pass = max_dist < error_tol;