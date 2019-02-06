function aug_array = augment_array_ndim(array,n)

dim = size(array,1);

for i  = 1:dim
    aug_array(i,:) = reshape([array(i,:);zeros(n-1,5)],[],1);
end

aug_array = kron(aug_array,ones(n,1));

for i = 1:dim
    for k = 1:n
       aug_array(n*(i-1) + k, :) = circshift(aug_array(n*(i-1) + k, :),k - 1); 
    end
end 