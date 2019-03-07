function [A_coll, b_coll, k_ctr] = ondemand_softconstraints(hor_rob,Phi,X0,A0,i,...
                                                 rmin,order,E1,E2)
A_coll = [];
b_coll = [];
K = size(hor_rob,2);
relax_lim = 10;
k_ctr = 0;
% Detect first collision on the horizon and return constraint tuple
for k = 1:K
    hor_k = squeeze(hor_rob(:,k,:));
    [collisions,neighbrs] = check_horizon(hor_k,i,E1,rmin,order); 
    if (any(collisions))
        k_ctr = k;
        N_v = sum(neighbrs);
        [A_coll,b_coll, dist] = build_constraint(hor_k,k,neighbrs,X0,Phi,...
                                                 A0,i,rmin,order,E1,E2);
        ncols = size(A_coll,2);
%         fprintf("k_ctr = %i with dist = %.2f m\n", k_ctr, dist)
        A_coll = [A_coll diag(dist);
                  zeros(N_v, ncols) eye(N_v);
                  zeros(N_v, ncols) -eye(N_v)];
%         if k_ctr <=6
%             relax_lim = 0.1;
%         end
        b_coll = [b_coll; zeros(N_v,1); relax_lim*ones(N_v,1)];
        return
    end
end