clc
clear all
close all
warning('off','all')
visualize = 0;  % 3D visualization of trajectory and predictions

% Define model parameters for the quad + controller system
model_params.zeta_xy = 0.6502;
model_params.tau_xy = 0.3815;
model_params.omega_xy = 1/model_params.tau_xy;
model_params.zeta_z = 0.9103;
model_params.tau_z = 0.3;
model_params.omega_z = 1/model_params.tau_z;

% Disturbance applied to the model within a time frame
disturbance = 1;
agent_disturb = [1];
disturbance_k = 40:60;

% dimension of space - 3 = 3D, 2 = 2D
ndim = 3; 

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
deg_poly = 3; % degree of differentiability required for the position
l = 3;  % number of Bezier curves to concatenate
d = 5;  % degree of the bezier curve

N = 1; % number of vehicles

% Noise standard deviation information based on Vicon data
std_p = 0.00228682;
std_v = 0.0109302;

% Physical limits
phys_limits.pmin = [-10,-10,0.2];
phys_limits.pmax = [10,10,4.2];
phys_limits.amax = 2;
phys_limits.amin = -2;

% Minimum distance between vehicles in m
rmin_init = 0.75;

% Initial positions
% [po,pf] = randomTest(N,pmin,pmax,rmin_init);

% Initial positions
po1 = [-2.0,2.0,1.5];
po2 = [-2.0,-2.01,1.5];
po3 = [-2.0,2.0,1.5];
po4 = [2.0,-2.0,1.5];
po = cat(3,po1,po2,po3,po4);

% Final positions
pf1 = [2.0,2.0,1.5];
pf2 = [2.0,2.0,1.5];
pf3 = [2.0,-2.0,1.5];
pf4 = [-2.0,2.0,1.5];
pf  = cat(3,pf1,pf2,pf3,pf4);

%% CONSTRUCT DOUBLE INTEGRATOR MODEL AND ASSOCIATED MATRICES
[model, inv_model] = get_model(h, model_params);
Lambda = get_lambda(model.A, model.B, k_hor);
A0 = get_a0(model.A, k_hor);

%% CONSTRUCT MATRICES TO WORK WITH BEZIER CURVES
% Beta - converts control points into polynomial coefficients
Beta = mat_bernstein2power(d,l,ndim);

% Gamma - sample the polynomial at different time steps
Gamma = mat_sample_poly(T_segment,0:h:((k_hor-1)*h),d,l);

% Alpha - sum of squared derivatives of the position
cr = zeros(1,d+1); % weights on degree of derivative deg = [0 1 ... d]
cr(3) = .008;
Alpha = mat_sum_sqrd_derivatives(d,T_segment,cr,l,ndim);
%% CONSTRUCT TERMS OF THE COST FUNCTION
% The Hessian for the minimum snap cost function is
H_snap = Beta'*Alpha*Beta;

% For the goal tracking error cost function, define a weight matrix S
s = 50;
spd = 3;
S = s*[zeros(3*(k_hor-spd), 3*k_hor);
       zeros(3*spd, 3*(k_hor-spd)) eye(3*spd)];
Phi = Lambda.pos*Gamma*Beta;
Phi_vel = Lambda.vel*Gamma*Beta;
H_err = Phi'*S*Phi;

% The complete Hessian is simply the sum of the two
H = H_err + H_snap;

% The linear term of the cost function depends both on the goal location of
% agent i and on its current position and velocity
% We can construct the static part that depends on the desired location
for i = 1:N
    f_pf(:,:,i) = repmat((pf(:,:,i))',k_hor,1)'*S*Phi;
end

% We can also construct the matrix that will then be multiplied by the
% initial condition -> X0'*A0'*S*Lambda*Gamma*Beta
mat_f_x0 = A0.pos'*S*Lambda.pos*Gamma*Beta;
%% CONSTRUCT INEQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Acceleration limits 2) workspace boundaries
% We need to create the matrices that map position control points into c-th
% derivative control points
T_ctrl_pts = cell_derivate_ctrl_pts(d);
[A_in, b_in] = build_ineq_constraints(d,l,h,ndim,k_hor,T_segment,...
                                      phys_limits,T_ctrl_pts, Beta);

%% CONSTRUCT EQUALITY CONSTRAINT MATRICES
% Types of constraints: 1) Continuity constraints up to degree deg_poly
A_eq = build_eq_constraints(d,l,ndim,deg_poly,T_ctrl_pts);
% The constant vector b_eq will be updated within the solver function,
% since it requires knowledge of the initial condition of the reference

%% MATRICES TO DECODE SOLUTION
% First construct all the matrices that map the solution vector to samples
% of the n-th derivative of position
[model_s, inv_model_s] = get_model(Ts,model_params);
K_sample = length(0:Ts:h-Ts);
Lambda_s = get_lambda(model_s.A,model_s.B,K_sample);
A0_s = get_a0(model_s.A,K_sample);

for r = 0:d
    if r > 0
        Mu = T_ctrl_pts{r};
        Mu_3d = augment_array_ndim(Mu,3);
        Sigma_r = kron(eye(l),Mu_3d);
    else
        Sigma_r = eye(3*(d+1)*l);
    end
    
    Beta_r = mat_bernstein2power(d-r,l,ndim);     
    
    % Sample Bezier curves at 1/h Hz for the whole horizon
    t_sample_r = 0:h:((k_hor-1)*h);
    Tau_r = mat_sample_poly(T_segment,t_sample_r,d-r,l);
    Der_h{r+1} = Tau_r*Beta_r*Sigma_r;
    
    % Sample Bezier curves at 1/Ts Hz for the first applied input
    t_sample_r = 0:Ts:h-Ts;
    Tau_r = mat_sample_poly(T_segment,t_sample_r,d-r,l);
    Der_ts{r+1} = Tau_r*Beta_r*Sigma_r;
end

% Sample states at 1/Ts Hz for the first applied input
t_sample_r = 0:Ts:h-Ts;
Tau_r = mat_sample_poly(T_segment,t_sample_r,d,l);
Phi_sample = Lambda_s.pos*Tau_r*Beta;
Phi_vel_sample = Lambda_s.vel*Tau_r*Beta;

%% INIT ALGORITHM
for i = 1:N
   poi = po(:,:,i)';
   voi = zeros(3,1);
   X0(:,i) = [poi; voi];
   pos_k_i(:,1,i) = poi;
   pos_k_i_sample(:,1,i) = poi;
   X0_ref(:,:,i) = [poi, voi, zeros(3,1) zeros(3,1) zeros(3,1)];
   prev_state(:,i) = X0(:,i);
   prev_input(:,i) = X0_ref(:,1,i);
   for r = 1:deg_poly+1
      ref(:,1,r,i) = X0_ref(:,r,i); 
      ref_sample(:,1,r,i) = X0_ref(:,r,i);
   end
   hor_ref(:,:,i,1) = repmat(poi,1,k_hor);
   hor_rob(:,:,i,1) = repmat(poi,1,k_hor);
end
pred_X0 = X0;

% Variables for reference replanning based on state feedback
err_tol_pos = 0.05;  % tolerance between predicted and sensed state
err_tol_vel = 0.5;
max_cost_value = 1;  % tolerance before resetting reference to state
ki = 0.0;  % integral term
integ_err(:,1) = zeros(3,1);
integ_err_vel(:,1) = zeros(3,1);


%% MAIN LOOP
for k = 2:K
    for i = 1:N 
        % Compare the expected and sensed position at time k
        err_pos(:,k) = X0(1:3,i) - pred_X0(1:3,i);
        err_vel(:,k) = X0(4:6,i) - pred_X0(4:6,i);
        
        % Compare the current position and the reference
        err_pos_ref(:,k) = X0(1:3,i) - X0_ref(:,1,i);
        err_vel_ref(:,k) = X0(4:6,i) - X0_ref(:,2,i);
        
        der_err_ref(:,k) = (err_pos_ref(:,k) - err_pos_ref(:,k-1))/h;
        
        cost(:,k) = (err_pos_ref(:,k).^5)./(-X0(4:6,i)+0.01);
        
        % Integral term on position
        integ_err(:,k) = integ_err(:,k-1) + err_pos_ref(:,k)*h;
        
        for n = 1:ndim
            % Reset integral term if close to ref
%             if abs(err_pos_ref(n,k)) < 0.5
%                 integ_err(n,k) = 0;
%             end
            
            % Filter noise on position for feedback
            if abs(err_pos(n,k)) < err_tol_pos
                err_pos(n,k) = 0;
            end
        end
        
        % Reset reference to state if the error grows large
        if any(abs(cost(:,k)) > max_cost_value*ones(3,1)) || any(cost(:,k) < -0.0005*ones(3,1))
            X0_ref(:,1,i) = X0(1:3,i);
            X0_ref(:,2,i) = X0(4:6,i);
            X0_ref(:,3:5,i) = zeros(3,3);
        else
            X0_ref(:,1,i) = X0_ref(:,1,i) + err_pos(:,k) + ki*integ_err(:,k);
        end
              
        f_tot = f_pf(:,:,i);

        % Include hard on demand collision avoidance
        [A_coll, b_coll] = ondemand_constraints(hor_rob(:,:,:,k-1),Phi,...
                                                X0(:,i),A0.pos,i,rmin,...
                                                order,E1,E2);
        % Augment the inequality constraints
        A_in_i = [A_in; A_coll];
        b_in_i = [b_in; b_coll];
        
        % Solve QP
        [x,exitflag] = MPC_update(l,deg_poly, A_in_i, b_in_i, A_eq, H,...
                                  mat_f_x0,f_tot, X0(:,i), X0_ref(:,:,i));
        if isempty(x)
            fprintf("ERROR: No solution - exitflag %i\n",exitflag)
            break;
        end
     
        % Apply input to model starting form our previous init condition
        pos_i = vec2mat(Phi*x + A0.pos*X0(:,i),3)';
        vel_i = vec2mat(Phi_vel*x + A0.vel*X0(:,i),3)';
        
        % Sample at a higher frequency the interval 0:Ts:h-Ts
        % This tells us what should be the value of our state after
        % sending the optimal commands if the model was perfect
        pos_i_sample = vec2mat(Phi_sample*x + A0_s.pos*X0(:,i),3)';
        vel_i_sample = vec2mat(Phi_vel_sample*x + A0_s.vel*X0(:,i),3)';
        
        % Sample the resulting reference Bezier curves at 1/h and 1/Ts
        % Get the next input to be applied 'X0_ref'
        cols = 2+(k-2)*(h/Ts):1+(k-1)*(h/Ts);
        for r = 1:d+1
            rth_ref(:,:,r) = vec2mat(Der_h{r}*x,3)';
            rth_ref_sample(:,:,r) = vec2mat(Der_ts{r}*x,3)';
            X0_ref(:,r,i) = rth_ref(:,2,r);
            ref(:,k,r,i) = rth_ref(:,2,r);
            ref_sample(:,cols,r,i) = rth_ref_sample(:,:,r);
        end
        % Simulate sending trajectories every Ts and applying at each time
        % step noise to the measurements and propagating the state forward
        X0_ex(:,1) = X0(:,i);
        for k_ex = 2:length(t_sample_r) + 1
            X0_ex(:, k_ex -1) = X0_ex(:, k_ex -1) + rnd_noise(std_p,std_v);
            X0_ex(:,k_ex) = model_s.A*X0_ex(:, k_ex-1) + ...
                            model_s.B*rth_ref_sample(:, k_ex-1, 1);
        end
        
        if ~disturbance || ~ismember(k,disturbance_k)
            % Initial conditions for next MPC cycle - based on sensing
            X0(:,i) = X0_ex(:, end);

            % Update agent's states at 1/h and 1/Ts frequencies
            pos_k_i_sample(:,cols,i) = X0_ex(1:3, 2:end);
            vel_k_i_sample(:,cols,i) = X0_ex(4:6, 2:end);
        elseif disturbance && ismember(k,disturbance_k)
            % APPLY SIMULATED DISTURBANCE
            X0(1,i) = X0(1,i)-0.1;
            X0(4,i) = -0.5;
            prev_state(:,i) = X0(:,i);
            pos_k_i_sample(:,cols,i) = repmat(X0(1:3,i),1,h/Ts);
            vel_k_i_sample(:,cols,i) = repmat(X0(4:6,i),1,h/Ts);
        end
        
        pred_X0(:,i) = [pos_i_sample(:,end); vel_i_sample(:,end)];
        pos_k_i(:,k,i) = X0(1:3,i);
        vel_k_i(:,k,i) = X0(4:6,i);            
        prev_input(:,i) = rth_ref_sample(:,end-1,1);
        
        % Reference and state prediction horizons - visualization purposes
        hor_ref(:,:,i,k) = rth_ref(:,:,1);
        hor_rob(:,:,i,k) = pos_i;
%         figure(1)
%         plot(hor_ref(1,:,i,k))
%         
    end
    if isempty(x)
       break;
    end
end
%% 3D VISUALIZATION
if visualize
    figure(1)
    colors = distinguishable_colors(N);
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'currentchar',' ')
    while get(gcf,'currentchar')==' '
        for i = 1:N
        h_line(i) = animatedline('LineWidth',2,'Color',colors(i,:),'LineStyle',':');
        end
        for k = 1:K
            for i = 1:N
                clearpoints(h_line(i));
                addpoints(h_line(i),hor_rob(1,:,i,k),hor_rob(2,:,i,k),hor_rob(3,:,i,k));     
                hold on;
                grid on;
                xlim([phys_limits.pmin(1), phys_limits.pmax(1)])
                ylim([phys_limits.pmin(2), phys_limits.pmax(2)])
                zlim([0,phys_limits.pmax(3)])
                plot3(pos_k_i(1,k,i),pos_k_i(2,k,i),pos_k_i(3,k,i),'o',...
                    'LineWidth',2,'Color',colors(i,:));
                plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                      'LineWidth',2,'Color',colors(i,:));
                plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
                      'LineWidth',2,'Color',colors(i,:));    
            end
        drawnow
        end
        clf
        pause(0.5)
    end
end

%% PLOT INTER-AGENT DISTANCES OVER TIME
% figure(6)
% for i = 1:N
%     for j = 1:N
%         if(i~=j)
%             differ = E1*(pos_k_i_sample(:,:,i) - pos_k_i_sample(:,:,j));
%             dist = (sum(differ.^order,1)).^(1/order);
%             plot(t, dist, 'LineWidth',1.5);
%             grid on;
%             hold on;
%             xlabel('t [s]')
%             ylabel('Inter-agent distance [m]');
%         end
%     end
% end
% plot(t,rmin*ones(length(t),1),'--r','LineWidth',1.5);
%% PLOT STATES AND REFERENCE TRAJECTORIES
state_label = {'x', 'y', 'z'};
der_label = {'p', 'v', 'a', 'j', 's'};

figure(1)
state = 1;
grid on;
hold on;
plot(t, pos_k_i_sample(state,:,1),'Linewidth',1.5)
plot(t, ref_sample(state,:,1,1),'--r','Linewidth',1.5)
ylabel([state_label{state} ' [m]'])
xlabel ('t [s]')

figure(2)
state = 1;
derivative = 2;
grid on;
hold on;
plot(t, vel_k_i_sample(state,:,1),'Linewidth',1.5)
plot(t, ref_sample(state,:,derivative,1),'--r','Linewidth',1.5)
ylabel([ der_label{derivative} state_label{state}  ' [m]'])
xlabel ('t [s]')

% figure(2)
% state = 2;
% grid on
% hold on;
% plot(tk, pos_k_i(state,:,1),'Linewidth',1.5)
% plot(tk, ref(state,:,1,1),'--r','Linewidth',1.5)
% ylabel([state_label{state} ' [m]'])
% xlabel ('t [s]')
% 
% figure(3)
% state = 3;
% grid on
% hold on;
% plot(tk, pos_k_i(state,:,1),'Linewidth',1.5)
% plot(tk, ref(state,:,1,1),'--r','Linewidth',1.5)
% ylabel([state_label{state} ' [m]'])
% xlabel ('t [s]')
% 

figure(3)
state = 1;
derivative = 3;
grid on
hold on;
plot(t, ref_sample(state,:,derivative,1),'Linewidth',1.5)
ylabel([ der_label{derivative} state_label{state}  ' [m]'])
xlabel ('t [s]')

figure(4)
subplot(3,1,1)
plot(tk, err_pos_ref(1,:),'Linewidth',1.5)
grid on
ylabel('Error [m]')
xlabel ('tk [s]')

subplot(3,1,2)
plot(tk, der_err_ref(1,:),'Linewidth',1.5)
grid on
ylabel('Derivative of Error [m]')
xlabel ('tk [s]')

subplot(3,1,3)
plot(tk, integ_err(1,:),'Linewidth',1.5)
grid on
ylabel('Integral of Error [m]')
xlabel ('tk [s]')

figure(5)
plot(tk, cost(1,:),'Linewidth',1.5)
grid on
ylabel('Cost')
xlabel ('tk [s]')



% figure(5)
% state = 1;
% derivative = 4;
% grid on
% hold on;
% plot(tk, ref(state,:,derivative,1),'Linewidth',1.5)
% ylabel([ der_label{derivative} state_label{state}  ' [m]'])
% xlabel ('t [s]')    