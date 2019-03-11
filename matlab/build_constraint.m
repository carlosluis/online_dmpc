function [Ain_k, bin_k, prev_dist, pf_tmp] = build_constraint(hor_k, k, neighbours, X0, Phi,...
                                                              A0, i, rmin, order, E1, E2)

N = size(hor_k, 2);
K = size(A0,1) / 3;
N_coll = sum(neighbours);
Ain_k = zeros(N_coll, size(Phi,2));
bin_k = zeros(N_coll, 1);
k_ctr = k-1;
p_i = X0(1:3);
idx = 1;
pf_tmp = [];
if k_ctr == 0
    k_ctr = 1;
end

% Loop through all the colliding neighbours and append the contraint tuple (A, b)
for j = 1:N
   if (i~= j && neighbours(j))
       p_j = hor_k(:,j);
       dist = norm(E1*(p_i-p_j), order);
       differ = (E2*(p_i-p_j).^(order-1))'; 
       prev_dist(idx) = dist^(order-1);
       
       if k_ctr <=1 && dist < rmin-0.05
           % We're almost crashing, force agents to repel each other
           pf_tmp = p_i + (p_i-p_j)*(rmin+0.1 -dist)/dist;
       end
       
       % Build intermediate terms for matrices
       init_cond = differ * A0(3 *(k_ctr-1) + 1 : 3*k_ctr, :)*X0;
       zeta = dist^(order-1)*(rmin - dist);
       rho = zeta + differ*p_i - init_cond;
       diff_mat = [zeros(1, 3*(k_ctr-1)) differ zeros(1, 3*(K-k_ctr))];
       
       % Build inequality constraint Ain_k*x <= bin_k
       Ain_k(idx,:) = -diff_mat*Phi;
       bin_k(idx) = -rho;
       idx = idx + 1;
   end
end
