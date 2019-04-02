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
N_vector = [5 10 15 20 25 30];
trials = 50;

for i = 1:N_vector(end)
    order(i) = order_a;
    rmin(i) = rmin_a;
    c(i,:) = c_a;
    E1(:,:,i) = E1_a;
    E2(:,:,i) = E2_a;
end

build_all_matrices;

pmin_gen = [-1.5,-1.5, 0.2];
pmax_gen = [1.5, 1.5, 2.2];

% Start Test
for p = 1:length(N_vector)
    N = N_vector(p);
    
    for  q = 1:trials
        fprintf("Doing trial #%i with %i vehicles\n", q, N)
        
        % Generate a random set of initial and final positions
        [po, pf] = random_test(N, pmin_gen, pmax_gen, rmin + 0.2, E1, order);
        
        % Run algorithm with on demand collision avoidance with repel strategy
        use_ondemand = true;
        use_repel = true;
        run_algorithm;
        
        % record the metrics for comparison
        repel_usage{p, q} = num_repels;
        avgrepel_usage(p, q) = mean(num_repels);
        dmpc_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
        dmpc_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
        dmpc_violated(p,q) = violated;
        dmpc_reachedgoal(p, q) = pass;
        dmpc_success(p, q) = ~violated && pass;
        if pass
            dmpc_totdist(p, q) = sum(sum(sqrt(diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 )));
                                          
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
        
         % Run algorithm with on demand collision avoidance without repel strategy
        use_ondemand = true;
        use_repel = false;
        run_algorithm;
        
        % record the metrics for comparison
        dmpc2_tbuild(p,q) = 1000*mean2(t_build(3:end,:));
        dmpc2_tqp(p,q) = 1000*mean2(t_qp(3:end,:));
        dmpc2_violated(p,q) = violated;
        dmpc2_reachedgoal(p, q) = pass;
        dmpc2_success(p, q) = ~violated && pass;
        if pass
            dmpc2_totdist(p, q) = sum(sum(sqrt(diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 )));
                                          
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
        
        % Run algorithm with BVC collision avoidance
        use_ondemand = false;
        run_algorithm;
        
        % record the metrics for comparison
        bvc_tbuild(p, q) = 1000*mean2(t_build(3:end,:));
        bvc_tqp(p, q) = 1000*mean2(t_qp(3:end,:));
        bvc_violated(p, q) = violated;
        bvc_reachedgoal(p, q) = pass;
        bvc_success(p, q) = ~violated && pass;
        if pass
            bvc_totdist(p, q) = sum(sum(sqrt(diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 + ...
                                              diff(ref_sample(1,:,1,:)).^2 )));
                                          
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
    end    
end
fprintf("Finished! \n")
save('comp_repel_bvc2')

%% Post Processing

% Probability of success plots
prob_dmpc = sum(dmpc_success, 2) / trials;
prob_dmpc2 = sum(dmpc2_success, 2) / trials;
prob_bvc = sum(bvc_success, 2) / trials;

figure(1)
grid on
hold on
plot(N_vector, prob_dmpc, 'r', 'Linewidth', 2)
plot(N_vector, prob_dmpc2, 'g', 'Linewidth', 2)
plot(N_vector, prob_bvc, 'b', 'Linewidth', 2)
xlabel('Number of agents');
ylabel('Success Probability');
legend('On-demand w/ repel', 'On-demand w/o repel','BVC')

% Plot times
figure(2)
grid on
hold on
tbuild_dmpc = mean(dmpc_tbuild, 2);
tqp_dmpc = mean(dmpc_tqp, 2);
tbuild_bvc = mean(bvc_tbuild, 2);
tqp_bvc = mean(bvc_tqp, 2);
plot(N_vector, tqp_dmpc, 'r', 'Linewidth', 2);
plot(N_vector, tbuild_dmpc, '--r', 'Linewidth', 2);
plot(N_vector, tqp_bvc, 'b', 'Linewidth', 2);
plot(N_vector, tbuild_bvc, '--b', 'Linewidth', 2);
xlabel('Number of agents');
ylabel('Computation Time per Agent [ms]');
legend('On-demand solve QP', 'On-demand build constr', 'BVC solve QP', 'BVC build constr')

% Completion time
tmean_traj = nanmean(dmpc_trajtime,2);
tstd_traj = nanstd(dmpc_trajtime,1,2);

tmean_traj2 = nanmean(dmpc2_trajtime,2);
tstd_traj2 = nanstd(dmpc2_trajtime,1,2);

tmean_traj3 = nanmean(bvc_trajtime,2);
tstd_traj3 = nanstd(bvc_trajtime,1,2);
figure(3)
grid on;
hold on;
errorbar(N_vector, tmean_traj, tstd_traj, 'r', 'Linewidth', 2);
errorbar(N_vector, tmean_traj2, tstd_traj2, 'g', 'Linewidth', 2);
errorbar(N_vector, tmean_traj3, tstd_traj3, 'b', 'Linewidth',2);
xlabel('Number of agents');
ylabel('Average Time for Transition [s]');
legend('On-demand w/ repel', 'On-demand w/o repel','BVC')

