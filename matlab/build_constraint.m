function [Ain_k, bin_k, prev_dist] = build_constraint(hor_k,k, neighbours, X0, Phi,...
                                                      A0,i,rmin,order,E1,E2)

N = size(hor_k,2);
K = size(A0,1)/3;
N_coll = sum(neighbours);
Ain_k = zeros(N_coll, size(Phi,2));
bin_k = zeros(N_coll,1);
k_ctr = k;
p_i = hor_k(:,i);
idx = 1;
if k_ctr == 0
    k_ctr = 1;
end

for j = 1:N
   if (i~= j && neighbours(j))
       p_j = hor_k(:,j);
       dist = norm(E1*(p_i-p_j),order);
       differ = (E2*(p_i-p_j).^(order-1))'; % Transpose of the difference
       
       prev_dist(idx) = dist^(order-1);
       
       % Build intermediate terms for matrices
       init_cond = differ*A0(3*(k_ctr-1)+1:3*k_ctr,:)*X0;
       zeta = dist^(order-1)*(rmin - dist);
       rho = zeta + differ*p_i - init_cond;
       diff_mat = [zeros(1,3*(k_ctr-1)) differ zeros(1,3*(K-k_ctr))];
       
       % Build inequality constraint Ain_k*x <= bin_k
       Ain_k(idx,:) = -diff_mat*Phi;
       bin_k(idx) = -rho;
    
       idx = idx + 1;
   end
end
