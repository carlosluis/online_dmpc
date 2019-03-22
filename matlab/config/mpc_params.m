% Clear variables
clear;

% Create a .mat file with the MPC tuning parameters

%%%%%%%%%% Energy minimization cost %%%%%%%%%%%%%%%%%%
cost_acc = .008;

%%%%%%%%%% Goal tracking error cost %%%%%%%%%%%%%%%%%% 
% s refers to cost and spd refers to the amount of time steps 
% to include in the error minimization
% Case 1: no collisions in the horizon - go fast
s_free = 100;
spd_f = 5;

% Case 2: collisions in the horizon - go slower
s_obs = 100;
spd_o = 1;

% Case 3: colliding scenario, change setpoint and repel rapidly
s_repel = 1000;
spd_r = 10;

%%%%%%%%%% Collision relaxation penalty %%%%%%%%%%%%
lin_coll_penalty = -1*10^8;
quad_coll_penalty = 1*10^6;

%%%%%%%%%% Tolerances to trigger a replan %%%%%%%%%%%%
err_tol_pos = 0.05;         % tolerance between predicted and sensed pos
err_tol_vel = 0.5;          % tolerance between predicted and sensed vel
max_cost = 0.8*ones(3, 1);   % max threshold on the strain cost function
min_cost = -0.01*ones(3, 1); % min threshold on the strain cost function
ki = 0.0;                   % integral term, currently not used

save('mpc_params.mat')

