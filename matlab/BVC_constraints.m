function [A_coll, b_coll] = BVC_constraints(current_pos, d, i,rmin,order,E1,E2,x_length)

% Construct the hyperspace constraints that will limit the first
% segment of the Bezier curve to lie within the BVC.

N = size(current_pos,2);
p_i = current_pos(:,i);
idx = 1;
A_coll = zeros(d*(N-1),x_length);
b_coll = zeros(d*(N-1),1);
for j = 1:N
   if (i~= j)
       p_j = current_pos(:,j);
       dist = norm(E1*(p_j-p_i),order);
       differ = (E2*(p_i-p_j).^(order-1))'; % Transpose of the difference
       
       % Right side of inequality constraint
       r = dist^(order-1)*((rmin/2+dist/2) - dist + differ*p_i/(dist^(order-1)));
       
       % the diff 1x3 vector must multiply all the control points of the
       % first segment of the Bezier curve, except the first one
       for  k = 1:d
          A_coll((idx-1)*d + k, 3*(k-1)+4:3*(k-1)+6) = -differ;
          b_coll((idx-1)*d + k) = -r;
       end
       idx = idx + 1;
   end
    
end
