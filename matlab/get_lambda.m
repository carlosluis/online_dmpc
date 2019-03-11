function Lambda = get_lambda(A, B, K)
Lambda.pos = [];
Lambda.vel = [];
prev_row = zeros(6, 3*K); 
ncols = size(B, 2);
nrows = size(B, 1);
for k = 1:K
    add_B = [zeros(nrows, ncols*(k-1)), B, zeros(nrows, ncols*(K-k))];
    new_row = A*prev_row + add_B;   
    Lambda.pos = [Lambda.pos; new_row(1:3,:)];
    Lambda.vel = [Lambda.vel; new_row(4:6,:)];
    prev_row = new_row;   
end