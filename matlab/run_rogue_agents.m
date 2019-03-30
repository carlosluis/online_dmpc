clc
clear all
close all
warning('off','all')

% Load simulation parameters
load('sim_params.mat')

% Load MPC tuning parameters
load('mpc_params.mat')

% Choose what data to visualize
visualize = 0;      % 3D visualization of trajectory and predictions
view_states = 1;    % pos, vel and acc of all agents
view_distance = 0;  % inter-agent distance over time
view_cost = 0;      % value of the replanning cost function
global debug_constr;
debug_constr = 0;

use_ondemand = true;

% Disturbance applied to the model within a time frame
disturbance = 1;       % activate the disturbance
agent_disturb = [1];   % choose which agents to perturb
disturbance_k = [1:20, 50:70];  % timesteps to apply the perturbation

% We will assume that all the rogue agents are labelled after the commanded agents

% Number of vehicles in the problem
N = 1;
N_rogues = 0;

% Specify a specific size for rogue agents
order_r = 2;
rmin_r = 0.75;
c_r = [1.0, 1.0, 2.0];
E_r = diag(c_r);
E1_r = E_r^(-1);
E2_r = E_r^(-order_r);

% Number of agents to be controlled by our algorithm
N_cmd = N - N_rogues;

for i = 1:N
    if i <= N_cmd
        order(i) = order_a;
        rmin(i) = rmin_a;
        c(i,:) = c_a;
        E1(:,:,i) = E1_a;
        E2(:,:,i) = E2_a;
    else
        order(i) = order_r;
        rmin(i) = rmin_r;
        c(i,:) = c_r;
        E1(:,:,i) = E1_r;
        E2(:,:,i) = E2_r;
        
    end
end

pmin_gen = [-1.5,-1.5,0.2];
pmax_gen = [1.5,1.5,2.2];

% Generate a random set of initial and final positions
% [po, pf] = random_test_static_rogues(N, N_cmd, pmin_gen, pmax_gen, rmin, E1, order);

% Initial positions
po1 = [0.0, 1.0,1.0];
po2 = [-1.0,-1.0,1.0];
po3 = [-1.0,1.0,1.0];
po4 = [1.0,-1.0,1.0];
po5 = [1.0, 0.0, 1.0];
po6 = [-1.0, 0.0, 1.0];
po7 = [-0.0, 0.0, 1.0];
po = cat(3,po1,po2,po7,po4,po5,po6,po7);
% 
% % Final positions
pf1 = [1.0, 1.0,1.0];
pf2 = [1.0,1.0,1.0];
pf3 = [1.0,-1.0,1.0];
pf4 = [-1.0,1.0,1.0];
pf5 = [-1.0, 0.0, 1.0];
pf6 = [1.0, 0.0, 1.0];
pf  = cat(3,pf1);

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

% The linear term of the cost function depends both on the goal location of
% agent i and on its current position and velocity
% We can construct the static part that depends on the desired location
for i = 1:N_cmd
    f_pf_free(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_free*Phi;
    f_pf_obs(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_obs*Phi;
end

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

%% %%%%%%%%%%%%%% INIT ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
   poi = po(:,:,i)';
   voi = 0.001*ones(3, 1);
   X0(:,i) = [poi; voi];
   pos_k_i(:,1,i) = poi;
   vel_k_i(:,1,i) = voi;
   pos_k_i_sample(:,1,i) = poi;
   X0_ref(:,:,i) = [poi, voi, zeros(3,d-1)];
   prev_state(:,i) = X0(:,i);
   for r = 1:deg_poly+1
      ref(:,1,r,i) = X0_ref(:,r,i); 
      ref_sample(:,1,r,i) = X0_ref(:,r,i);
   end
   hor_ref(:,:,i,1) = repmat(poi, 1, k_hor);
   hor_rob(:,:,i,1) = repmat(poi, 1, k_hor+1);
end

pred_X0 = X0;

% Variables for reference replanning based on state feedback
integ_err(:,1) = zeros(3, 1);
colors = distinguishable_colors(N);
tic
%% %%%%%%%%%%%%% MPC MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_build = zeros(K, N_cmd);
% t_qp = zeros(K, N_cmd);
for k = 2:K
    
    % Update states for the regular agents
    for i = 1:N_cmd 
        % Compare the expected and sensed position at time k
        err_pos(:,k) = X0(1:3,i) - pred_X0(1:3,i);
        err_vel(:,k) = X0(4:6,i) - pred_X0(4:6,i);
        
        % Compare the current position and the reference
        err_pos_ref(:,k) = X0(1:3,i) - X0_ref(:,1,i);
        err_vel_ref(:,k) = X0(4:6,i) - X0_ref(:,2,i);
        
        der_err_ref(:,k) = (err_pos_ref(:,k) - err_pos_ref(:,k-1)) / h;
        
        % Cost that determines whether there's something disturbing the agent
        % Cost gets higher when the error between the reference and the state gets higher
        cost(:,k,i) = (err_pos_ref(:,k).^5) ./ (-X0(4:6,i)+sign(X0(4:6,i))*0.01);
        
        % Integral term on position
        integ_err(:,k) = integ_err(:,k-1) + err_pos_ref(:,k)*h;
        
        for n = 1:ndim
            % Filter noise on position for feedback
            if abs(err_pos(n,k)) < err_tol_pos
                err_pos(n,k) = 0;
            end
        end
        
        trigger(k,i) = 0;
        
        % Reset reference to state if the error grows large
        if any(cost(:,k,i) > max_cost) || any(cost(:,k,i) < min_cost)
            X0_ref(:,1,i) = X0(1:3,i);
            X0_ref(:,2,i) = X0(4:6,i);
            X0_ref(:,3:5,i) = zeros(3,3);
            trigger(k,i) = 1;
        else
            X0_ref(:,1,i) = X0_ref(:,1,i); %+ err_pos(:,k) + ki*integ_err(:,k);
        end
              
        % Include on-demand collision avoidance
        
        if use_ondemand
            [A_coll, b_coll, pf_tmp, t_build(k,i)] = ondemand_softconstraints(hor_rob(:,2:end,:,k-1), Phi,...
                                                            X0(:,i), A0.pos, i, rmin,...
                                                            order, E1, E2);

            if ~isempty(b_coll) % collisions in the horizon
                % Include collision constraints and slack variables
                N_v = length(b_coll) / 3;
                A_in_i = [A_in zeros(size(A_in,1), N_v) ; A_coll];
                b_in_i = [b_in; b_coll];
                A_eq_i = [A_eq zeros(size(A_eq,1), N_v)];

                % Linear and quadratic term to penalize collision relaxation
                f_eps = lin_coll_penalty*ones(1, N_v);
                H_eps = quad_coll_penalty*eye(N_v);

                % If close to colliding, change setpoint to quickly react
                if ~isempty(pf_tmp)
                    H_i = [H_r zeros(size(H_f,1), N_v);
                           zeros(N_v, size(H_f,2)) H_eps];
                    mat_f_x0_i = mat_f_x0_repel;
                    f_tot = repmat((pf_tmp),k_hor,1)'*Rho_repel;
                else
                    H_i = [H_o zeros(size(H_f,1), N_v);
                           zeros(N_v,size(H_f,2)) H_eps];
                    mat_f_x0_i = mat_f_x0_obs;
                    f_tot = f_pf_obs(:,:,i);
                end

            else % no collisions in horizon
                A_in_i = A_in;
                b_in_i = b_in;
                A_eq_i = A_eq;
                H_i = H_f;
                f_eps = [];
                mat_f_x0_i = mat_f_x0_free;
                f_tot = f_pf_free(:,:,i);
            end
        else % Use BVC constraints
            x_length = (d+1) * ndim * l;
            t_start = tic;
            [A_coll, b_coll] = BVC_constraints_ref(X0_ref, d, i, rmin, order, E1, E2, x_length);
            t_build(k,i) = toc(t_start);
            A_in_i = [A_in; A_coll];
            b_in_i = [b_in; b_coll];
            A_eq_i = A_eq;
            H_i = H_f;
            f_eps = [];
            mat_f_x0_i = mat_f_x0_free;
            f_tot = f_pf_free(:,:,i);
        end
        
        % Solve QP
        t_start = tic;
        [sol, exitflag] = softMPC_update(l, deg_poly, A_in_i, b_in_i, A_eq_i, H_i,...
                                      mat_f_x0_i, f_tot, f_eps, X0(:,i), X0_ref(:,:,i));
                                  
        t_qp(k,i) = toc(t_start);  
        if  isempty(sol)
            fprintf("No solution found, using previous input \n")
            x = prev_x{i};
%             assert(~isempty(x), 'ERROR: No solution found - exitflag =  %i\n',exitflag);
        else
            prev_x{i} = sol;
            x = sol;
        end
        
        % Extract the control points
        u = x(1:size(mat_f_x0_free, 2));
     
        % Apply input to model starting form our previous init condition
        pos_i = vec2mat(Phi*u + A0.pos*X0(:,i),3)';
        vel_i = vec2mat(Phi_vel*u + A0.vel*X0(:,i),3)';
        
        if ~isempty(b_coll) && debug_constr && use_ondemand
            figure(1)
            plot3(pos_i(1,:), pos_i(2,:), pos_i(3,:),...
                  '*','Color',colors(i,:),'Linewidth',2)
            plot3(pf(1,1,i),pf(1,2,i),pf(1,3,i),'s','Color',colors(i,:),'Linewidth',4,'markers',10)
            fprintf("relaxed constraint by %.2f cm\n", epsilon*100)
            hola = 1;
        end
        
        % Sample at a higher frequency the interval 0:Ts:h-Ts
        % This tells us what should be the value of our state after
        % sending the optimal commands if the model was perfect
        pos_i_sample = vec2mat(Phi_sample*u + A0_s.pos*X0(:,i),3)';
        vel_i_sample = vec2mat(Phi_vel_sample*u + A0_s.vel*X0(:,i),3)';
        
        % Sample the resulting reference Bezier curves at 1/h and 1/Ts
        % Get the next input to be applied 'X0_ref'
        cols = 2 + (k-2)*(h/Ts):1 + (k-1)*(h/Ts);
        for r = 1:d+1
            rth_ref(:,:,r) = vec2mat(Der_h{r}*u, 3)';
            rth_ref_sample(:,:,r) = vec2mat(Der_ts{r}*u, 3)';
            X0_ref(:,r,i) = rth_ref(:,2,r);
            ref(:,k,r,i) = rth_ref(:,2,r);
            ref_sample(:,cols,r,i) = rth_ref_sample(:,:,r);
        end
        
        % Simulate sending trajectories every Ts and applying at each time
        % step noise to the measurements and propagating the state forward
        X0_ex(:,1) = X0(:,i);
        for k_ex = 2:length(t_sample_r) + 1
            X0_ex(:, k_ex -1) = X0_ex(:, k_ex -1) + rnd_noise(std_p,std_v);
            X0_ex(:,k_ex) = model_s.A*X0_ex(:, k_ex-1) + model_s.B*rth_ref_sample(:, k_ex-1, 1);
        end
        
        if ~disturbance || ~ismember(k,disturbance_k) || ~ismember(i,agent_disturb)
            % Initial conditions for next MPC cycle - based on sensing
            X0(:,i) = X0_ex(:, end);
            
            % Update agent's states at 1/h and 1/Ts frequencies
            pos_k_i_sample(:,cols,i) = X0_ex(1:3, 2:end);
            vel_k_i_sample(:,cols,i) = X0_ex(4:6, 2:end);
            
        elseif disturbance && ismember(k,disturbance_k) && ismember(i,agent_disturb)
            % Apply simulated disturbance
            if k <= 30
                X0(1,i) = X0(1,i);
                X0(4,i) = 0;
            elseif k >= 50
                X0(1,i) = X0(1,i) - 0.05;
                X0(4,i) = -0.25;
            end
                
            X0(:,i) = X0(:,i) + rnd_noise(std_p,std_v);
            prev_state(:,i) = X0(:,i);
            pos_k_i_sample(:,cols,i) = repmat(X0(1:3,i),1,h/Ts);
            vel_k_i_sample(:,cols,i) = repmat(X0(4:6,i),1,h/Ts);
        end
        
        pred_X0(:,i) = [pos_i_sample(:,end); vel_i_sample(:,end)];
        pos_k_i(:,k,i) = X0(1:3,i);
        vel_k_i(:,k,i) = X0(4:6,i);            
        
        % Reference and state prediction horizons - visualization purposes
        hor_ref(:,:,i,k) = rth_ref(:,:,1);
        hor_rob_k(:,:,i) = [X0(1:3,i) pos_i(:,1:end)];
        
        if i==0
            hor_rob_k(:,:,i) = [X0(1:3,i) repmat(X0(1:3,i),1,k_hor)];
        else
            hor_rob_k(:,:,i) = [X0(1:3,i) pos_i(:,1:end)];
        end
    end
  
    % Update the states for the rogue agents
    for  i = N_cmd + 1:N
        % For now we're forcing the agents to remain in place
        pos_k_i(:,k,i) = pos_k_i(:,k-1,i);
        vel_k_i(:,k,i) = vel_k_i(:,k-1,i);
        hor_rob_k(:,:,i) = repmat(X0(1:3, i), 1, k_hor + 1);
    end

    hor_rob(:,:,:,k) = hor_rob_k;
end
toc

%% POST-CHECK OF SOLUTION

% Check if collision constraints were not violated
violated = false;
rmin_check = 0.15;
c_check = 3;
E_check = diag([1,1,c_check]);
E1_check = E_check^(-1);

for i = 1:N_cmd
    for j = 1:N
        if(i~=j)
            differ = E1_check*(pos_k_i(:,:,i) - pos_k_i(:,:,j));
            dist = (sum(differ.^order(j),1)).^(1/order(j));
            if min(dist) < rmin_check
                violated = true;
                [value,index] = min(dist);
                fprintf("Collision violation by %.2fcm: vehicles %i and %i @ t = %.2fs \n",...
                        (rmin(j) -value)*100,i,j,index*h);
            end
        end
    end
end

if ~violated
    fprintf("No collisions during execution!\n");
end
    
% Check if all vehicles reached their goals.
pass = reached_goal(pos_k_i(:,:,1:N_cmd), pf, 0.1, N_cmd);
if pass
    fprintf("All agents reached their goals\n");
else
    fprintf("Some agents didn't reached their goals\n");
end

% Print average building time and qp solving time for each agent
fprintf("Average collision constraint building time = %.2f ms\n", 1000*mean2(t_build(3:end,:)));
fprintf("Average QP solving time = %.2f ms\n", 1000*mean2(t_qp(3:end,:)));


%% %%%%%%%%%%% PLOT STATES AND REFERENCE TRAJECTORIES %%%%%%%%%%%%%%%%%%%%%%%%%%%

state_label = {'x', 'y', 'z'};
der_label = {'p', 'v', 'a', 'j', 's'};
% colors = distinguishable_colors(N);

if view_states
    for i = 1:N_cmd
        % Position
        figure(1)
        subplot(3,1,1)
        state = 1;
        grid on;
        hold on;
        plot(t, pos_k_i_sample(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, ref_sample(state,:,1,i), '--', 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, phys_limits.pmin(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.pmax(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel([state_label{state} ' [m]'])
        xlabel ('t [s]')

        subplot(3,1,2)
        state = 2;
        grid on;
        hold on;
        plot(t, pos_k_i_sample(state,:,i), 'Linewidth',1.5, 'Color', colors(i,:))
        plot(t, ref_sample(state,:,1,i), '--', 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, phys_limits.pmin(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.pmax(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel([state_label{state} ' [m]'])
        xlabel ('t [s]')

        subplot(3,1,3)
        state = 3;
        grid on;
        hold on;
        plot(t, pos_k_i_sample(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, ref_sample(state,:,1,i), '--', 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, phys_limits.pmin(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.pmax(state)*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel([state_label{state} ' [m]'])
        xlabel ('t [s]')
        
        % Velocity
        figure(2)
        subplot(3,1,1)
        state = 1;
        grid on;
        hold on;
        plot(t, vel_k_i_sample(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        ylabel(['v' state_label{state} ' [m/s]'])
        xlabel ('t [s]')

        subplot(3,1,2)
        state = 2;
        grid on;
        hold on;
        plot(t, vel_k_i_sample(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        ylabel(['v' state_label{state} ' [m/s]'])
        xlabel ('t [s]')

        subplot(3,1,3)
        state = 3;
        grid on;
        hold on;
        plot(t, vel_k_i_sample(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        ylabel(['v' state_label{state} ' [m/s]'])
        xlabel ('t [s]')
        
        % Acceleration of the reference
        figure(3)
        subplot(3,1,1)
        state = 1;
        grid on;
        hold on;
        plot(t, ref_sample(state,:,3,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, phys_limits.amin*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.amax*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel(['a' state_label{state} ' [m/s^2]'])
        xlabel ('t [s]')

        subplot(3,1,2)
        state = 2;
        grid on;
        hold on;
        plot(t, ref_sample(state,:,3,i),'Linewidth',1.5, 'Color', colors(i,:))
        plot(t, phys_limits.amin*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.amax*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel(['a' state_label{state} ' [m/s^2]'])
        xlabel ('t [s]')

        subplot(3,1,3)
        state = 3;
        grid on;
        hold on;
        plot(t, ref_sample(state,:,3,i), 'Linewidth', 1.5, 'Color', colors(i,:))
        plot(t, phys_limits.amin*ones(length(t),1), '--r', 'LineWidth', 1.5);
        plot(t, phys_limits.amax*ones(length(t),1), '--r', 'LineWidth', 1.5);
        ylabel(['a' state_label{state} ' [m/s^2]'])
        xlabel ('t [s]')
        
        % Replanning threshold cost function
        if view_cost
            figure(4)
            subplot(3,1,1)
            state = 1;
            grid on;
            hold on;
            plot(tk, cost(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
            ylabel(['Cost in ' state_label{state}])
            xlabel ('t [s]')

            subplot(3,1,2)
            state = 2;
            grid on;
            hold on;
            plot(tk, cost(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
            ylabel(['Cost in ' state_label{state}])
            xlabel ('t [s]')

            subplot(3,1,3)
            state = 3;
            grid on;
            hold on;
            plot(tk, cost(state,:,i), 'Linewidth', 1.5, 'Color', colors(i,:))
            ylabel(['Cost in ' state_label{state}])
            xlabel ('t [s]')
        end
    end
end

%% %%%%%%%%%% PLOT INTER-AGENT DISTANCES OVER TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if view_distance
    figure(6)
    for i = 1:N
        for j = 1:N
            if(i~=j)
                differ = E1*(pos_k_i(:,:,i) - pos_k_i(:,:,j));
                dist = (sum(differ.^order,1)).^(1/order);
                plot(tk, dist, 'LineWidth', 1.5);
                grid on;
                hold on;
                xlabel('t [s]')
                ylabel('Inter-agent distance [m]');
            end
        end
    end
    plot(tk, rmin*ones(length(tk), 1), '--r', 'LineWidth', 1.5);
end

%% %%%%%%%%%%%% 3D VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xi, Yi, Zi] = sphere;
if visualize
    figure(1)
    colors = distinguishable_colors(N);
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'currentchar',' ')
    while get(gcf,'currentchar')==' '
        for i = 1:N
        h_line(i) = animatedline('LineWidth', 2, 'Color', colors(i,:), 'LineStyle', ':');
        end
        for k = 1:K
            for i = 1:N_cmd
                if k ~= 1
                    delete(h_pos(i))
                end
                clearpoints(h_line(i));
                addpoints(h_line(i), hor_rob(1,:,i,k), hor_rob(2,:,i,k), hor_rob(3,:,i,k));
                hold on;
                grid on;
                xlim([phys_limits.pmin(1), phys_limits.pmax(1)])
                ylim([phys_limits.pmin(2), phys_limits.pmax(2)])
                zlim([0, phys_limits.pmax(3)+1.0])
                h_pos(i) = plot3(pos_k_i(1,k,i), pos_k_i(2,k,i), pos_k_i(3,k,i), 'o',...
                                 'LineWidth', 2, 'Color',colors(i,:));
                plot3(po(1,1,i), po(1,2,i), po(1,3,i), '^',...
                      'LineWidth', 2, 'Color', colors(i,:));
                plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i), 'x',...
                      'LineWidth', 2, 'Color', colors(i,:));   
            end
            for i = N_cmd + 1: N
                % Plot rouge agents sphere for better visualization
                XX = Xi * rmin(i) + pos_k_i(1,k,i);
                YY = Yi * rmin(i) + pos_k_i(2,k,i);
                ZZ = Zi * rmin(i) * c(i) + pos_k_i(3,k,i);
                surface(XX, YY, ZZ, 'Facealpha', 0.5, 'FaceColor',...
                        [0.3,0.3,0.3], 'EdgeColor', [0,0,0]);
                
%                 plot3(pos_k_i(1,k,i), pos_k_i(2,k,i), pos_k_i(3,k,i), 'k*',...
%                                  'LineWidth', 2, 'Markers', 12);
            end
            
            
        drawnow
        end
        pause(0.5)
        clf
    end
end