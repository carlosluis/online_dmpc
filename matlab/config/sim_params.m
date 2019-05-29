% Create a .mat file with all fixed simulation parameters to solve
% the optimization, as well as the tuning parameters

% Clear everything
clear;

% Model parameters for the quad + controller system
% It's model as a discrete second order system
model_params.zeta_xy = 0.6502;
model_params.tau_xy = 0.3815;
model_params.omega_xy = 1/model_params.tau_xy;
model_params.zeta_z = 0.9103;
model_params.tau_z = 0.3;
model_params.omega_z = 1/model_params.tau_z;

% VICON Measurements noise std for position and velocity data
std_p = 0.00228682;
std_v = 0.0109302;

% Dimension of space - 3 = 3D, 2 = 2D
ndim = 3; 

% Time settings and variables
T = 20;          % simulation duration
h = 0.2;         % time step between MPC updates
tk = 0:h:T;      % coarse discrete time base
K = T/h + 1;     % number of time steps to simulate
Ts = 0.01;       % send h/Ts cmds in-between MPC updates 
t = 0:Ts:T;      % interpolated time vector
k_hor = 16;      % horizon length - duration of (k_hor-1)*h sec
T_segment = 1.0; % fixed time length of each Bezier segment

% Collision ellipsoid parameters
order_a = 2;         % order of the ellipsoid - choose between 2 and 4
rmin_a = 0.3;       % X-Y protection radius for collisions
c_a = [1.0, 1.0, 2.0];           % Z radius is equal to rmin*c
E_a = diag(c_a); % scaling vector to compute distances to ellipsoid 
E1_a = E_a^(-1);
E2_a = E_a^(-order_a);

% Bezier curve parameters. Note that d > deg_poly always
deg_poly = 3;  % degree of differentiability required for the position
l = 3;         % number of Bezier curves to concatenate
d = 5;         % degree of the bezier curve

% Physical limits of the robot - position and acceleration bounds
phys_limits.pmin = [-1.5, -1.5, 0.98];
phys_limits.pmax = [1.5, 1.5, 1.02];
phys_limits.amax = 1;
phys_limits.amin = -1;

save('sim_params.mat')
