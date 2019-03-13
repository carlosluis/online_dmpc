function aug_array = augment_array_ndim(array, n)

nrows = size(array, 1);
ncols = size(array, 2);

for i  = 1:nrows
    aug_array(i,:) = reshape([array(i,:); zeros(n-1,ncols)], [], 1);
end

aug_array = kron(aug_array, ones(n,1));

for i = 1:nrows
    for k = 1:n
       aug_array(n*(i-1) + k, :) = circshift(aug_array(n*(i-1) + k, :), k-1); 
    end
end 
hola = 1;