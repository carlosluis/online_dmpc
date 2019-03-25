function [A_coll, b_coll] = BVC_constraints_state(X0, Phi, A0, i, rmin, order, E1, E2,...
                                                  x_length, K, k_coll)

% Construct the hyperspace constraints that will limit the first
% segment of the Bezier curve to lie within the BVC.
N = size(X0,2);
p_i = X0(1:3,i);
idx = 1;
A_coll = zeros(k_coll*(N-1), x_length);
b_coll = zeros(k_coll*(N-1), 1);
for j = 1:N
   if (i~= j)
       p_j = X0(1:3,j);
       dist = norm(E1*(p_j-p_i),order);
       differ = (E2*(p_i-p_j).^(order-1))'; % Transpose of the difference
       
       for k = 1:k_coll
           % RHS of inequality constraint
           r = dist^(order-1)*((rmin/2+dist/2) - dist + differ*p_i/(dist^(order-1))) - differ*A0(3*(k-1)+1:3*k,:)*X0;
           diff_mat = [zeros(1,3*(k-1)) differ zeros(1,3*(K-k))];
           % LHS of inequality constraint Ain*x <= bin
           A_coll((idx-1)*k_coll + k,:) = -diff_mat*Phi;
           b_coll((idx-1)*k_coll + k) = -r;
       end
       idx = idx + 1;
   end
end
