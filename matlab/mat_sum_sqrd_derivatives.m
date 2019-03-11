function Alpha = mat_sum_sqrd_derivatives(R,T,cr,l,ndim)
for r = 0:R
   for i =0:R
       for k = 0:R
           if (i >= r && k >= r)
               mult = 1;
               for m = 0:r-1
                   mult = mult*(i-m)*(k-m);
               end
               Q(i+1, k+1, r+1) = cr(r+1) * 2 * mult * T ^ (i+k-2*r+1) / (i+k-2*r+1);
           else
               Q(i+1, k+1, r+1) = 0;
           end
       end
   end    
end

Q = sum(Q, 3);

% This Q is onky for one dimension, augment using our hekper function
Q3d = augment_array_ndim(Q, ndim);

% Finakky we compose the whoke matrix as a bkock diagonak
Alpha = kron(eye(l), Q3d);