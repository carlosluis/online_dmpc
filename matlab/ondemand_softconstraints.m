function [A_coll, b_coll, pf_tmp] = ondemand_softconstraints(hor_rob, Phi, X0, A0, i,...
                                                             rmin, order, E1, E2)
A_coll = [];
b_coll = [];
K = size(hor_rob, 2);
relax_lim = 2;
pf_tmp = [];
% Detect first collision on the horizon and return constraint tuple
for k = 1:K
    hor_k = squeeze(hor_rob(:,k,:));
    [collisions, neighbrs] = check_horizon(hor_k, i, E1, rmin, order); 
    if (any(collisions))
        k_ctr = k;
        N_v = sum(neighbrs);
        [A_coll, b_coll, dist, pf_tmp] = build_constraint(hor_k, k, neighbrs, X0, Phi,...
                                                          A0, i, rmin, order, E1, E2);
        ncols = size(A_coll, 2);
%         fprintf("k_ctr = %i with dist = %.2f m \n", k_ctr, min(dist));
        A_coll = [A_coll diag(dist);
                  zeros(N_v, ncols) eye(N_v);
                  zeros(N_v, ncols) -eye(N_v)];

        b_coll = [b_coll; zeros(N_v, 1); relax_lim*ones(N_v, 1)];
        return
    end
end