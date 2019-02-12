clc
clear all
close all
warning('off','all')

% Time settings and variables
T = 20; % Trajectory final time
h = 0.2; % time step duration
tk = 0:h:T;
K = T/h + 1; % number of time steps
Ts = 0.01; % period for interpolation @ 100Hz
t = 0:Ts:T; % interpolated time vector
k_hor = 16; % horizon length
T_segment = 1.0; % fixed time length of each Bezier segment

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.5; % X-Y protection radius for collisions
c = 1.5; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Bezier curve parameters. Note that d > deg_poly always
deg_poly = 2; % degree of differentiability required for the position
l = 3;  % number of Bezier curves to concatenate
d = 5;  % degree of the bezier curve

N = 1; % number of vehicles

% Workspace boundaries
pmin = [-2.5,-2.5,0.2];
pmax = [2.5,2.5,2.2];

% Acceleration limits
amax = 1;
amin = -1;

% Minimum distance between vehicles in m
rmin_init = 0.75;

% Initial positions
% [po,pf] = randomTest(N,pmin,pmax,rmin_init);

% Initial positions
po1 = [1.0,1.0,1.5];
po2 = [-1.5,-1.5,1.5];
po3 = [-1.5,1.5,1.5];
po4 = [1.5,-1.5,1.5];
po = cat(3,po1);

% Final positions
pf1 = [2.0,1.0,1.5];
pf2 = [1.5,1.5,1.5];
pf3 = [1.5,-1.5,1.5];
pf4 = [-1.5,1.5,1.5];
pf  = cat(3,pf1);

%% CONSTRUCT DOUBLE INTEGRATOR MODEL AND ASSOCIATED MATRICES

% Define model parameters for the quad + controller system
zeta_xy = 0.6502;
tau_xy = 0.3815;
omega_xy = 1/tau_xy;
zeta_z = 0.9103;
tau_z = 0.3;
omega_z = 1/tau_z;

% Matrices of the discrete state space model
% State = [x y z vx vy vz]
A = [1  0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     -h*omega_xy^2 0 0 -2*omega_xy*h*zeta_xy 0 0;
     0 -h*omega_xy^2 0 0 -2*omega_xy*h*zeta_xy 0;
     0 0 -h*omega_xy^2 0 0 -2*omega_xy*h*zeta_xy];
 
B = [zeros(3,3);
     h*omega_xy^2 0 0;
     0 h*omega_xy^2 0;
     0 0 h*omega_xy^2]; 
 
Lambda = getLambda(A,B,k_hor);
Lambda_K = Lambda(end-2:end,:); % to be used in the cost function

A0 = getA0(A,k_hor);
A0_K = A0(end-2:end,:);

%% CONSTRUCT MATRICES TO WORK WITH BEZIER CURVES
% Delta - converts control points into polynomial coefficients
delta = Bern2power(d);
delta3d = augment_array_ndim(delta,3);

% Then we assemble a block-diagonal matrix based on the number of segments
% this matrix is used to compute the minimum snap cost function
Beta = kron(eye(l),delta3d);

% Beta is also used to sample the Bezier curve at different times
% We will use it to minimize the error at the end of the curve, by sampling
% our curve at discrete time steps with length 'h'
Tau = getTauPos(h,T_segment,d);
Tau3d_all = augment_array_ndim(Tau,3);

% To eliminate overlapping control points between segments
Tau3d_k = Tau3d_all(1:end-3,:);

Gamma = blkdiag(kron(eye(l-1),Tau3d_k),Tau3d_all);

% Construct matrix Q: Hessian for each of the polynomial derivatives
cr = zeros(1,d+1); % weights on degree of derivative deg = [0 1 ... d]
cr(5) = .01;
Q = getQ(d,T_segment,cr);
Q = sum(Q,3);

% This Q is only for one dimension, augment using our helper function
Q3d = augment_array_ndim(Q,3);

% Finally we compose the whole matrix as a block diagonal
Alpha = kron(eye(l),Q3d);

%% CONSTRUCT TERMS OF THE COST FUNCTION
% The Hessian for the minimum snap cost function is
H_snap = Beta'*Alpha*Beta;

% For the trajectory tracking error cost function, define a weight matrix S
s = 100;
spd = 1;
S = s*[zeros(3*(k_hor-spd),3*k_hor);
       zeros(3*spd,3*(k_hor-spd)) eye(3*spd)];
Phi = Lambda*Gamma*Beta;
H_err = Phi'*S*Phi;

% The complete Hessian is simply the sum of the two
H = H_snap + H_err;

% The linear term of the cost function depends both on the goal location of
% agent i and on its current position and velocity
% We can construct the static part that depends on the desired location
for i = 1:N
    f_pf(:,:,i) = repmat((pf(:,:,i))',k_hor,1)'*S*Lambda*Gamma*Beta;
end

% We can also construct the matrix that will then be multiplied by the
% initial condition -> X0'*A0'*S*Lambda*Gamma*Beta
mat_f_x0 = A0'*S*Lambda*Gamma*Beta;
%% CONSTRUCT INEQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Acceleration limits 2) workspace boundaries

% We need to create the matrices that map position control points into c-th
% derivative control points
der_mats = getDerivatives(d);
n = d; % n represents the n-th derivative of position
for k = 1:d
    aux = 1;
    for k_aux = d:-1:k
       aux = der_mats{k_aux}*aux;   
    end
    T_ctrl_pts{n}(:,:) = aux; 
    n = n - 1;
end

% ACCELERATION CONSTRAINT
% Acceleration is a quadratic Bezier curve, The mapping between position
% and acceleration control points is given by the linear mapping
M = T_ctrl_pts{2};

% convert this matrix into a 3d version
M_3d = augment_array_ndim(M,3);
Sigma_acc = kron(eye(l),M_3d);

% Multiplying the decision vector X times Sigma_acc gives acc control pts
% We define a delta_acc that converts acc control points into polynomials
delta_acc = Bern2power(d-2);
         
delta_acc_3d = augment_array_ndim(delta_acc,3);
Beta_acc = kron(eye(l),delta_acc_3d);

% Define Tau_acc that evaluates the polynomial at different times
% Define the times at which we desire to sample between 0s and 2.4s
t_sample_acc = 0:(h):((k_hor-1)*h);
Tau_acc = getTauSamples(T_segment,t_sample_acc,d-2,l);

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_acc = [Tau_acc*Beta_acc*Sigma_acc;
           -Tau_acc*Beta_acc*Sigma_acc];

% The vector is defined by the maximum and minimum acc allowed 
b_in_acc = [amax*ones(3*length(t_sample_acc),1);
           -amin*ones(3*length(t_sample_acc),1)];

% POSITION REFERENCE CONSTRAINT - WORKSPACE BOUNDARIES
t_sample_pos_ref = 0:(h):((k_hor-1)*h);
Tau_pos_ref = getTauSamples(T_segment,t_sample_pos_ref,d,l);

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_pos_ref = [Tau_pos_ref*Beta;
                -Tau_pos_ref*Beta];

% The vector b_in_pos_ref depends solely on the limits pmin and pmax
b_in_pos_ref = [repmat(pmax',length(t_sample_pos_ref),1);
                repmat(-pmin',length(t_sample_pos_ref),1)];
    

% The complete matrix for constraints is just the concatenation of both
A_in = [A_in_acc; A_in_pos_ref];
b_in = [b_in_acc; b_in_pos_ref];

%% CONSTRUCT EQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Continuity constraints up to degree deg_poly
% From the control points of each derivative of the Bezier curve, we must
% select the first point and the overlapping points between segments
D = zeros(l, 3*(d+1));
for k = 1:l
    if k==1
       D(1,1) = 1; 
    else
       D(k,(k-1)*(d+1)) = 1;
       D(k,(k-1)*(d+1) + 1) = -1;
    end
end

A_eq = augment_array_ndim(D,3);

% Make T_ctrl_pts to be a 3D matrix and representing l Bezier segments 
if deg_poly > 0
    T_ctrl_pts_3d{deg_poly} = [];
    D_der{deg_poly} = [];
    for k = 1:deg_poly
        % Select the ctrl points we need and sustract them
        for n = 1:l
            if n == 1
                D_der{k}(1,:) = [T_ctrl_pts{k}(1,:) zeros(1,(l-1)*(d+1))];
            else
                cols = (n-2)*(d+1)+1: (n)*(d+1);
                D_der{k}(n,cols) = [T_ctrl_pts{k}(end,:) -T_ctrl_pts{k}(1,:)];
            end
        end
        % Construct equality constraint matrix by stacking matrices together
        A_eq = [A_eq;augment_array_ndim(D_der{k},3)];
    end
end
% The constant vector b_eq will be updated within the solver function,
% since it requires knowledge of the initial condition of the reference

%% MAIN LOOP
% We initialize the algorithm with every agent in its initial position and
% every derivative to be zero.

for k = 1:K
    for i = 1:N
        if k == 1
            poi = po(:,:,i)';
            voi = zeros(3,1);
            X0 = [poi; voi];
            X0_ref = [poi voi 0*ones(3,1) zeros(3,1) zeros(3,1)];
            x = MPC_update(l,deg_poly, A_in, b_in, A_eq, H, mat_f_x0, f_pf(:,:,i), X0, X0_ref);
            pos = Phi*x + A0*X0;
            t_sample_bla = 0:0.01:((k_hor-1)*h);
            Tau_bla = getTauSamples(T_segment,t_sample_bla,d,l);
            pos_ref = Tau_bla*Beta*x;
            p_ref = vec2mat(pos_ref,3)';
            p = vec2mat(pos,3)';
            
            figure(1)
            plot(0:(1*h):((k_hor-1)*h),p(1,:))
            hold on
            plot(t_sample_bla,p_ref(1,:))
            
        end      
    end
end

