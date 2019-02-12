function delta = Bern2power(n)
% Converts a polynomial in the bernstein basis to a polynomial in the
% regular power basis [1 t t^2 ... t^n]
delta = zeros(n+1,n+1);
for k = 0:n
   for i = k:n
       delta(i+1,k+1) = (-1)^(i-k)*nchoosek(n,i)*nchoosek(i,k);
   end
end