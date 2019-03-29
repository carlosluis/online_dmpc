%%%%%%%%%%%%%% CONSTRUCT DOUBLE INTEGRATOR MODEL AND ASSOCIATED MATRICES %%%%%%%%%

[model, inv_model] = get_model(h, model_params);
Lambda = get_lambda(model.A, model.B, k_hor);
A0 = get_a0(model.A, k_hor);

%%%%%%%%%%%%% CONSTRUCT MATRICES TO WORK WITH BEZIER CURVES %%%%%%%%%%%%%%%%%%%%%

% Beta - converts control points into polynomial coefficients
Beta = mat_bernstein2power(d, l, ndim);

% Gamma - sample the polynomial at different time steps
Gamma = mat_sample_poly(T_segment, 0:h:((k_hor-1)*h), d, l);

% Alpha - sum of squared derivatives of the position
cr = zeros(1, d+1); % weights on degree of derivative deg = [0 1 ... d]
cr(3) = cost_acc;
Alpha = mat_sum_sqrd_derivatives(d, T_segment, cr, l, ndim);

%%%%%%%%%%%%%% CONSTRUCT TERMS OF THE COST FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Hessian for the minimum snap cost function is
H_snap = Beta'*Alpha*Beta;

% For the goal tracking error cost function, define a weight matrix S
% Case 1: no collisions in the horizon - go fast
S_free = s_free*[zeros(3*(k_hor-spd_f), 3*k_hor);
                 zeros(3*spd_f, 3*(k_hor-spd_f)) eye(3*spd_f)];
             
% Case 2: collisions in the horizon - go slower
S_obs = s_obs*[zeros(3*(k_hor-spd_o), 3*k_hor);
               zeros(3*spd_o, 3*(k_hor-spd_o)) eye(3*spd_o)];

% Case 3: colliding scenario, change setpoint and repel rapidly
S_repel = s_repel*[eye(3*spd_r) zeros(3*spd_r, 3*(k_hor-spd_r));
                   zeros(3*(k_hor-spd_r), 3*k_hor)];
   
Phi = Lambda.pos*Gamma*Beta;
Phi_vel = Lambda.vel*Gamma*Beta;
H_free = Phi'*S_free*Phi;
H_obs = Phi'*S_obs*Phi;
H_repel = Phi'*S_repel*Phi;

% The complete Hessian is simply the sum of the two
H_f = H_free + H_snap;
H_o = H_obs + H_snap;
H_r = H_repel + H_snap;

% Predefine this matrix to construct f_pf_repel when needed
Rho_repel = S_repel*Phi;

% We can also construct the matrix that will then be multiplied by the
% initial condition -> X0'*A0'*S*Lambda*Gamma*Beta
mat_f_x0_free = A0.pos'*S_free*Lambda.pos*Gamma*Beta;
mat_f_x0_obs = A0.pos'*S_obs*Lambda.pos*Gamma*Beta;
mat_f_x0_repel = A0.pos'*S_repel*Lambda.pos*Gamma*Beta;

%%%%%%%%%%%%% CONSTRUCT INEQUALITY CONSTRAINT MATRICES %%%%%%%%%%%%%%%%%%%%%%

% Types of constraints: 1) Acceleration limits 2) workspace boundaries
% We need to create the matrices that map position control points into c-th
% derivative control points
T_ctrl_pts = cell_derivate_ctrl_pts(d);
[A_in, b_in] = build_ineq_constraints(d, l, h, ndim, k_hor, T_segment,...
                                      phys_limits, T_ctrl_pts, Beta);

%%%%%%%%%%%%% CONSTRUCT EQUALITY CONSTRAINT MATRICES %%%%%%%%%%%%%%%%%%%%

% Types of constraints: 1) Continuity constraints up to degree deg_poly
A_eq = build_eq_constraints(d, l, ndim, deg_poly, T_ctrl_pts);
% The constant vector b_eq will be updated within the solver function,
% since it requires knowledge of the initial condition of the reference

%%%%%%%%%%%%%% MATRICES TO DECODE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First construct all the matrices that map the solution vector to samples
% of the n-th derivative of position
[model_s, inv_model_s] = get_model(Ts, model_params);
K_sample = length(0:Ts:h-Ts);
Lambda_s = get_lambda(model_s.A, model_s.B, K_sample);
A0_s = get_a0(model_s.A, K_sample);

for r = 0:d
    if r > 0
        Mu = T_ctrl_pts{r};
        Mu_3d = augment_array_ndim(Mu, 3);
        Sigma_r = kron(eye(l), Mu_3d);
    else
        Sigma_r = eye(3*(d+1)*l);
    end
    
    Beta_r = mat_bernstein2power(d-r, l, ndim);     
    
    % Sample Bezier curves at 1/h Hz for the whole horizon
    t_sample_r = 0:h:((k_hor-1)*h);
    Tau_r = mat_sample_poly(T_segment, t_sample_r, d-r, l);
    Der_h{r+1} = Tau_r*Beta_r*Sigma_r;
    
    % Sample Bezier curves at 1/Ts Hz for the first applied input
    t_sample_r = 0:Ts:h-Ts;
    Tau_r = mat_sample_poly(T_segment, t_sample_r, d-r, l);
    Der_ts{r+1} = Tau_r*Beta_r*Sigma_r;
end

% Sample states at 1/Ts Hz for the first applied input
t_sample_r = 0:Ts:h-Ts;
Tau_r = mat_sample_poly(T_segment, t_sample_r, d, l);
Phi_sample = Lambda_s.pos*Tau_r*Beta;
Phi_vel_sample = Lambda_s.vel*Tau_r*Beta;
