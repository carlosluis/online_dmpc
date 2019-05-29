function [A_in,b_in] = build_ineq_constraints(d,l,h,ndim,K,T_segment,...
                                              phys_limits,T_ctrl_pts,Beta)

%%%%%%%%%%% ACCELERATION CONSTRAINT %%%%%%%%%%%%%%%%%%%%%%%%%

% Acceleration is a quadratic Bezier curve, The mapping between position
% and acceleration control points is given by the linear mapping
M = T_ctrl_pts{2};

% Convert this matrix into a 3d version
M_3d = augment_array_ndim(M, ndim);
Sigma_acc = kron(eye(l), M_3d);

% Multiplying the decision vector X times Sigma_acc gives acc control pts
% We define a delta_acc that converts acc control points into polynomials
Beta_acc = mat_bernstein2power(d-2, l, ndim);

% Define Tau_acc that evaluates the polynomial at different times
% Define the times at which we desire to sample between 0s and 2.4s
% NOTE: Always constrain the first 
t_sample_acc = h:(2*h): (K-1)*h;
Tau_acc = mat_sample_poly(T_segment, t_sample_acc, d-2, l);

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_acc = [Tau_acc*Beta_acc*Sigma_acc;
           -Tau_acc*Beta_acc*Sigma_acc];

% The vector is defined by the maximum and minimum acc allowed 
b_in_acc = [phys_limits.amax*ones(3*length(t_sample_acc), 1);
           -phys_limits.amin*ones(3*length(t_sample_acc), 1)];

% POSITION REFERENCE CONSTRAINT - WORKSPACE BOUNDARIES
t_sample_pos_ref = h:(h) : (K-1)*h;
Tau_pos_ref = mat_sample_poly(T_segment, t_sample_pos_ref, d, l);

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_pos_ref = [Tau_pos_ref*Beta;
                -Tau_pos_ref*Beta];

% The vector b_in_pos_ref depends solely on the limits pmin and pmax
b_in_pos_ref = [repmat(phys_limits.pmax', length(t_sample_pos_ref), 1);
                repmat(-phys_limits.pmin', length(t_sample_pos_ref), 1)];
    
% The complete matrix for constraints is just the concatenation of both
% A_in = [A_in_acc];
% b_in = [b_in_acc];

% A_in = [A_in_acc; A_in_pos_ref];
% b_in = [b_in_acc; b_in_pos_ref];

A_in = [A_in_acc; A_in_pos_ref];
b_in = [b_in_acc; b_in_pos_ref];