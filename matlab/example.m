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
k_hor = 13; % horizon length

% Variables for ellipsoid constraint
order = 2; % choose between 2 or 4 for the order of the super ellipsoid
rmin = 0.5; % X-Y protection radius for collisions
c = 1.5; % make this one for spherical constraint
E = diag([1,1,c]);
E1 = E^(-1);
E2 = E^(-order);

% Degree of differentiability required
deg_poly = 4;

N = 4; % number of vehicles

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
po1 = [1.501,1.5,1.5];
po2 = [-1.5,-1.5,1.5];
po3 = [-1.5,1.5,1.5];
po4 = [1.5,-1.5,1.5];
po = cat(3,po1,po2,po3,po4);

% Final positions
pf1 = [-1.5,-1.5,1.5];
pf2 = [1.5,1.5,1.5];
pf3 = [1.5,-1.5,1.5];
pf4 = [-1.5,1.5,1.5];
pf  = cat(3,pf1,pf2,pf3,pf4);

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
     0 0 -h*omega_z^2 0 0 -2*omega_z*h*zeta_z];
 
B = [zeros(3,3);
     h*omega_xy^2 0 0;
     0 h*omega_xy^2 0;
     0 0 h*omega_z^2]; 
 
Lambda = getLambda(A,B,k_hor);
Lambda_K = Lambda(end-2:end,:); % to be used in the cost function

A0 = getA0(A,k_hor);
A0_K = A0(end-2:end,:);

%% CONSTRUCT MATRICES TO WORK WITH BEZIER CURVES
% Delta - converts control points into polynomial coefficients
delta = [1 0 0 0 0;
         -4 4 0 0 0 ;
         6 -12 6 0 0;
         -4 12 -12 4 0;
         1 -4 6 -4 1 ];
     
delta3d = augment_array_ndim(delta,3);

% Then we assemble a block-diagonal matrix based on the number of segments
% this matrix is used to compute the minimum snap cost function
l = 3; 
Beta = kron(eye(l),delta3d);

% Beta is also used to sample the Bezier curve at different times
% We will use it to minimize the error at the end of the curve, by sampling
% our curve at discrete time steps with length 'h'
T_segment = 0.8; % fixed time length of each Bezier segment
Tau = getTauPos(h,T_segment,deg_poly);
Tau3d_all = augment_array_ndim(Tau,3);

% To eliminate overlapping control points between segments
Tau3d_k = Tau3d_all(1:end-3,:);

Gamma = blkdiag(kron(eye(l-1),Tau3d_k),Tau3d_all);

% Construct matrix Q: Hessian for each of the polynomial derivatives
cr = [0 1 0 0 0];
Q = getQ(4,T_segment,cr);
Q = sum(Q,3);

% This Q is only for one dimension, augment using our helper function
Q3d = augment_array_ndim(Q,3);

% Finally we compose the whole matrix as a block diagonal
Alpha = kron(eye(l),Q3d);

%% CONSTRUCT TERMS OF THE COST FUNCTION
% The Hessian for the minimum snap cost function is
H_snap = Beta'*Alpha*Beta;

% For the trajectory tracking error cost function, define a weight matrix S
s = 10;
spd = 1;
S = s*[zeros(3*(k_hor-spd),3*k_hor);
       zeros(3*spd,3*(k_hor-spd)) eye(3*spd)];
Phi = Lambda*Gamma*Beta;
H_err = Phi'*S*Phi;

% The complete Hessian is simply the sum of the two
H = H_snap + H_err;

% For the linear term of the cost function, we define
chi = S*Lambda*Gamma*Beta;

%% CONSTRUCT INEQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Acceleration limits 2) workspace boundaries

% ACCELERATION CONSTRAINT
% Acceleration is a quadratic Bezier curve, The mapping between position
% and acceleration control points is given by the linear mapping
M = 3*[4 -8  4  0 0;
       0  4 -8  4 0;
       0  0  4 -4 4];
% convert this matrix into a 3d version
M_3d = augment_array_ndim(M,3);
Sigma_acc = kron(eye(l),M_3d);

% Multiplying the decision vector X times Sigma_acc gives acc control pts
% We define a delta_acc that converts acc control points into polynomials
delta_acc = [1 0 0;
             -2 2 0;
             1 -2 1];
         
delta_acc_3d = augment_array_ndim(delta_acc,3);
Beta_acc = kron(eye(l),delta_acc_3d);

% Define Tau_acc that evaluates the polynomial at different times
% Define the times at which we desire to sample between 0s and 2.4s
t_sample_acc = 0:(4*h):((k_hor-1)*h);
Tau_acc = getTauSamples(T_segment,t_sample_acc,deg_poly-2,l);
Gamma_acc = getGammaSamples(Tau_acc,length(t_sample_acc));

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_acc = [Gamma_acc*Beta_acc*Sigma_acc;
           -Gamma_acc*Beta_acc*Sigma_acc];

% The vector is defined by the maximum and minimum acc allowed 
b_in_acc = [amax*ones(3*length(t_sample_acc),1);
           -amin*ones(3*length(t_sample_acc),1)];

% POSITION REFERENCE CONSTRAINT - WORKSPACE BOUNDARIES
t_sample_pos_ref = 0:(4*h):((k_hor-1)*h);
Tau_pos_ref = getTauSamples(T_segment,t_sample_pos_ref,deg_poly,l);
Gamma_pos_ref = getGammaSamples(Tau_pos_ref,length(t_sample_pos_ref));

% Now we compose the constraint in (A,b) form to pass to the solver
A_in_pos_ref = [Gamma_pos_ref*Beta;
                -Gamma_pos_ref*Beta];

% The vector b will be defined within the solver function since it depends
% on the goal location of the i-th agent

% The complete matrix for constraints is just the concatenation of both
A_in = [A_in_acc;A_in_pos_ref];

%% CONSTRUCT EQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Continuity constraints up to degree deg_poly


