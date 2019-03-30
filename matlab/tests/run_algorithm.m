%%
%%%%%%%%%%%%%%%% INIT ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:N
   f_pf_free(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_free*Phi;
   f_pf_obs(:,:,i) = repmat((pf(:,:,i))', k_hor, 1)'*S_obs*Phi; 
   poi = po(:,:,i)';
   voi = 0.001*ones(3, 1);
   X0(:,i) = [poi; voi];
   pos_k_i(:,1,i) = poi;
   vel_k_i(:,1,i) = voi;
   pos_k_i_sample(:,1,i) = poi;
   X0_ref(:,:,i) = [poi, voi, zeros(3,d - 1)];
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

%%%%%%%%%%%%%%% MPC MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_build = zeros(K, N_cmd);
% t_qp = zeros(K, N_cmd);
for k = 2:K
    
    % Update states for the regular agents
    for i = 1:N
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
        
        % Initial conditions for next MPC cycle - based on sensing
        X0(:,i) = X0_ex(:, end);

        % Update agent's states at 1/h and 1/Ts frequencies
        pos_k_i_sample(:,cols,i) = X0_ex(1:3, 2:end);
        vel_k_i_sample(:,cols,i) = X0_ex(4:6, 2:end);
        
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

    hor_rob(:,:,:,k) = hor_rob_k;
end


%% POST-CHECK OF SOLUTION

% Check if collision constraints were not violated
violated = false;
rmin_check = 0.15;
c_check = 3;
E_check = diag([1,1,c_check]);
E1_check = E_check^(-1);

for i = 1:N
    for j = 1:N
        if(i~=j)
            differ = E1_check*(pos_k_i(:,:,i) - pos_k_i(:,:,j));
            dist = (sum(differ.^order(j),1)).^(1/order(j));
            if min(dist) < (rmin_check)
                violated = true;
                break;
            end
        end
    end
    if violated
        break;
    end
end

% Check if all vehicles reached their goals.
pass = reached_goal(pos_k_i(:,:,1:N), pf, 0.2, N);

