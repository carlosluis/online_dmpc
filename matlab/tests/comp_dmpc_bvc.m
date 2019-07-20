clc
clear all
close all
warning('off','all')

% Load simulation parameters
load('sim_params.mat')

% Load MPC tuning parameters
load('mpc_params.mat')

% Choose what data to visualize
global debug_constr;
debug_constr = 0;

% Testing parameters
N_vector = 4:4:20;
trials = 50;

for i = 1:N_vector(end)
    order(i) = order_a;
    rmin(i) = rmin_a;
    c(i,:) = c_a;
    E1(:,:,i) = E1_a;
    E2(:,:,i) = E2_a;
end

build_all_matrices;

pmin_gen = [-0.7937,-0.7937,0.2];
pmax_gen = [0.7937,0.7937,1.7874];

% pmin_gen = [-1.5, -1.5, 0.2];
% pmax_gen = [1.5, 1.5, 2.2];

% Start Test
for p = 1:length(N_vector)
    N = N_vector(p);
    
    for  q = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n", q, N)
        
        % Generate a random set of initial and final positions
        [po, pf] = random_test(N, pmin_gen, pmax_gen, rmin, E1, order);
        
        % Run algorithm using hard BVC collision avoidance
        use_ondemand = false;
        use_repel = false;
        use_stateCA = false;
        use_softBVC = false;
        deg_poly = 2;
        run_algorithm;
        
        % record the metrics for comparison
        bvc_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
        bvc_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
        bvc_violated(p,q) = violated;
        bvc_reachedgoal(p, q) = pass;
        bvc_success(p, q) = ~violated && pass;
        if pass
            bvc_totdist(p, q) = sum(sum(sqrt(diff(pos_k_i_sample(1,:,:)).^2 + ...
                                              diff(pos_k_i_sample(2,:,:)).^2 + ...
                                              diff(pos_k_i_sample(3,:,:)).^2 )));
                                          
            for i = 1:N
                diff_goal = pos_k_i_sample(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                hola = find(dist_goal >= 0.1, 1, 'last');
                if isempty(hola)
                    time_index(i) = 0;
                else
                    time_index(i) = hola + 1;
                end
            end
            bvc_trajtime(p, q) = max(time_index)*Ts;
            
        else
            bvc_totdist(p, q) = nan;
            bvc_trajtime(p, q) = nan;
            
        end
        
        % Run algorithm using soft BVC collision avoidance
        use_ondemand = false;
        use_repel = false;
        use_stateCA = false;
        use_softBVC = true;
        deg_poly = 2;
        run_algorithm;
        
        % record the metrics for comparison
        bvc2_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
        bvc2_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
        bvc2_violated(p,q) = violated;
        bvc2_reachedgoal(p, q) = pass;
        bvc2_success(p, q) = ~violated && pass;
        if pass
            bvc2_totdist(p, q) = sum(sum(sqrt(diff(pos_k_i_sample(1,:,:)).^2 + ...
                                              diff(pos_k_i_sample(2,:,:)).^2 + ...
                                              diff(pos_k_i_sample(3,:,:)).^2 )));
                                          
            for i = 1:N
                diff_goal = pos_k_i_sample(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                hola = find(dist_goal >= 0.1, 1, 'last');
                if isempty(hola)
                    time_index(i) = 0;
                else
                    time_index(i) = hola + 1;
                end
            end
            bvc2_trajtime(p, q) = max(time_index)*Ts;
            
        else
            bvc2_totdist(p, q) = nan;
            bvc2_trajtime(p, q) = nan;
            
        end
        
        % Run algorithm with Collision avoidance in the state space.
        use_ondemand = true;
        use_repel = false;
        use_stateCA = false;
        use_softBVC = false;
        deg_poly = 2;
        run_algorithm;
        
        % record the metrics for comparison
        dmpc_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
        dmpc_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
        dmpc_violated(p,q) = violated;
        dmpc_reachedgoal(p, q) = pass;
        dmpc_success(p, q) = ~violated && pass;
        if pass
            dmpc_totdist(p, q) = sum(sum(sqrt(diff(pos_k_i_sample(1,:,:)).^2 + ...
                                              diff(pos_k_i_sample(2,:,:)).^2 + ...
                                              diff(pos_k_i_sample(3,:,:)).^2 )));
                                          
            for i = 1:N
                diff_goal = pos_k_i_sample(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                hola = find(dist_goal >= 0.1, 1, 'last');
                if isempty(hola)
                    time_index(i) = 0;
                else
                    time_index(i) = hola + 1;
                end
            end
            dmpc_trajtime(p, q) = max(time_index)*Ts;
            
        else
            dmpc_totdist(p, q) = nan;
            dmpc_trajtime(p, q) = nan;
            
        end
        
        % Run algorithm with Collision avoidance in the input space.
        use_ondemand = true;
        use_repel = false;
        use_stateCA = false;
        use_softBVC = false;
        deg_poly = 2;
        run_algorithm;
        
        % record the metrics for comparison
        dmpc2_tbuild(p, q) = 1000*mean2(t_build(3:end,:));
        dmpc2_tqp(p, q) = 1000*mean2(t_qp(3:end,:));
        dmpc2_violated(p, q) = violated;
        dmpc2_reachedgoal(p, q) = pass;
        dmpc2_success(p, q) = ~violated && pass;
        if pass
            dmpc2_totdist(p, q) = sum(sum(sqrt(diff(pos_k_i_sample(1,:,:)).^2 + ...
                                              diff(pos_k_i_sample(2,:,:)).^2 + ...
                                              diff(pos_k_i_sample(3,:,:)).^2 )));
                                          
            for i = 1:N
                diff_goal = pos_k_i_sample(:,:,i) - repmat(pf(:,:,i),length(t),1)';
                dist_goal = sqrt(sum(diff_goal.^2,1));
                hola = find(dist_goal >= 0.1, 1, 'last');
                if isempty(hola)
                    time_index(i) = 0;
                else
                    time_index(i) = hola + 1;
                end
            end
            dmpc2_trajtime(p, q) = max(time_index)*Ts;
            
        else
            dmpc2_totdist(p, q) = nan;
            dmpc2_trajtime(p, q) = nan;
        end
    end    
end
fprintf("Finished! \n")
save('comp_allCA')

%% Post Processing

% Probability of success plots
prob_bvc = sum(bvc_success, 2) / trials;
prob_bvc2= sum(bvc2_success, 2) / trials;
prob_dmpc = sum(dmpc_success, 2) / trials;
prob_dmpc2 = sum(dmpc2_success, 2) / trials;

figure(1)
grid on
hold on
plot(N_vector, prob_dmpc, 'r', 'Linewidth', 2)
plot(N_vector, prob_dmpc2, 'b', 'Linewidth', 2)
plot(N_vector, prob_bvc, 'm', 'Linewidth', 2)
plot(N_vector, prob_bvc2, 'g', 'Linewidth', 2)
xlabel('Number of agents');
ylabel('Success Probability');
legend('DMPC State CA','DMPC Input CA', 'BVC Hard', 'BVC Soft')

% Plot times
figure(2)
grid on
hold on
tqp_dmpc = mean(dmpc_tqp, 2);
tqp_dmpc2 = mean(dmpc2_tqp, 2);
tqp_bvc = mean(bvc_tqp, 2);
tqp_bvc2 = mean(bvc2_tqp, 2);
plot(N_vector, tqp_dmpc, 'r', 'Linewidth', 2);
plot(N_vector, tqp_dmpc2, 'b', 'Linewidth', 2);
plot(N_vector, tqp_bvc, 'm', 'Linewidth', 2);
plot(N_vector, tqp_bvc2, 'g', 'Linewidth', 2);
xlabel('Number of agents');
ylabel('Computation Time per Agent [ms]');
legend('DMPC State CA','DMPC Input CA', 'BVC Hard', 'BVC Soft')

% Completion time
tmean_traj = nanmean(dmpc_trajtime,2);
tstd_traj = nanstd(dmpc_trajtime,1,2);
tmean_traj2 = nanmean(dmpc2_trajtime,2);
tstd_traj2 = nanstd(dmpc2_trajtime,1,2);
tmeanbvc_traj = nanmean(bvc_trajtime,2);
tstdbvc_traj = nanstd(bvc_trajtime,1,2);
tmeanbvc_traj2 = nanmean(bvc2_trajtime,2);
tstdbvc_traj2 = nanstd(bvc2_trajtime,1,2);
figure(3)
grid on;
hold on;
errorbar(N_vector, tmean_traj, tstd_traj, 'r', 'Linewidth', 2);
errorbar(N_vector, tmean_traj2, tstd_traj2, 'b', 'Linewidth',2);
errorbar(N_vector, tmeanbvc_traj, tstdbvc_traj, 'm', 'Linewidth', 2);
errorbar(N_vector, tmeanbvc_traj2, tstdbvc_traj2, 'g', 'Linewidth',2);
xlabel('Number of agents');
ylabel('Average Time for Transition [s]');
legend('DMPC State CA','DMPC Input CA', 'BVC Hard', 'BVC Soft')