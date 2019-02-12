function [Lambda, Lambda_vel] = getLambda(A,B,K)

Lambda = []; % local inequality constraint variable
Lambda_vel = [];
prev_row = zeros(6,3*K); % For the first iteration of constructing matrix Ain

% Build matrix to convert acceleration to position
for k = 1:K
    add_B = [zeros(size(B,1),size(B,2)*(k-1)) B zeros(size(B,1),size(B,2)*(K-k))];
    new_row = A*prev_row + add_B;   
    Lambda = [Lambda; new_row(1:3,:)];
    Lambda_vel = [Lambda_vel; new_row(4:6,:)];
    prev_row = new_row;   
end