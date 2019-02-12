function Gamma = getGammaSamples(Tau,num_samples,dim)
l = length(Tau);
N = dim + 1;
Gamma = zeros(3*num_samples,3*N*l); %ncols = 3D * N ctrl pts * l
Tau_3d{l} = [];
curr_row = 1;
for i = 1:l
    if ~isempty(Tau{i})
        Tau_3d{i} = augment_array_ndim(Tau{i},3);
        nrow = size(Tau_3d{i},1);
        ncol = size(Tau_3d{i},2);
        Gamma(curr_row:curr_row+nrow-1,(i-1)*ncol+1:i*ncol) = Tau_3d{i};
        curr_row = curr_row+nrow;
    end
end