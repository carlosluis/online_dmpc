function x = MPC_update(l, deg_poly, Ain, bin, Aeq, H, mat_f_tot, f_pf, X0, X0_ref)
constr_tol = 1e-3;
options = optimoptions('quadprog','ConstraintTolerance',constr_tol);

% Construct linear term of the function based on X0 
% X0' = [px py pz vx vy vz]
f = -2*(f_pf + X0'*mat_f_tot);

% Construct the vector beq based on X0_ref
% X0_ref = [px vx ax jx sx]
%          |py vy ay jy sy|
%          [pz vz az jz sz]
beq = zeros(size(Aeq,1),1);
for i = 1:deg_poly+1
    beq((i-1)*l +1 : i*l) = X0_ref(:,i);
end

% Solve the QP
[x,fval,exitflag] = quadprog(2*H,f',[],[],Aeq,beq,[],[],[],options);
