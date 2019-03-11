function [A_coll, b_coll] = ondemand_constraints(hor_rob, Phi, X0, A0, i,...
                                                 rmin, order, E1, E2)
A_coll = [];
b_coll = [];
K = size(hor_rob, 2);
% Detect first collision on the horizon and return constraint tuple
for k = 1:K
    hor_k = squeeze(hor_rob(:,k,:));
    [collisions, neighbrs] = check_horizon(hor_k, i, E1, rmin, order);
    if (any(collisions))
        [A_coll,b_coll] = build_constraint(hor_k, k, neighbrs, X0, Phi, A0, i,...
                                           rmin, order, E1, E2);
        return;
    end
end