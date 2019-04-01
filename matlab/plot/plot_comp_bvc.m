clc
clear
close all

load('comp_bvc5.mat')

% Probability of success plots
prob_dmpc = sum(dmpc_success, 2) / trials * 100;
prob_bvc = sum(bvc_success, 2) / trials * 100;

figure(1)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([5:5:30]);
ylim([0,105])
xlim([3,32])
h1 = plot(N_vector, prob_bvc,':b', 'Linewidth',2.5);
plot(N_vector, prob_bvc,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);
h2 = plot(N_vector, prob_dmpc,':r','Linewidth',2.5);
plot(N_vector, prob_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);
xlabel('Number of Agents')
ylabel('Success Probability [%]');
[h, icons, plots, s] = legend([h1,h2],'BVC', 'On-demand');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');



% Completion time
tmean_traj = nanmean(dmpc_trajtime,2);
tstd_traj = nanstd(dmpc_trajtime,1,2);
tmean_traj2 = nanmean(bvc_trajtime,2);
tstd_traj2 = nanstd(bvc_trajtime,1,2);

figure(2)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([5:5:30]);
xlim([3,32])
ylim([0,20])
h1 = plot(N_vector, tmean_traj2,':b','LineWidth', 2.5);
plot(N_vector, tmean_traj2,'ob', 'MarkerFaceColor', 'b','Linewidth',1.5,'markers',10);
h2 = plot(N_vector, tmean_traj,':r','LineWidth', 2.5);
plot(N_vector, tmean_traj,'or', 'MarkerFaceColor', 'r','Linewidth',1.5,'markers',10);
xlabel('Number of Agents')
ylabel('Transition Time [s]');
[h, icons, plots, s] = legend([h1,h2],'BVC', 'On-demand');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');


% Plot times
tbuild_dmpc = mean(dmpc_tbuild, 2);
tqp_dmpc = mean(dmpc_tqp, 2);
tbuild_bvc = mean(bvc_tbuild, 2);
tqp_bvc = mean(bvc_tqp, 2);
figure(3)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([5:5:30]);
xlim([3,32])
h1 = plot(N_vector, tqp_bvc, ':b', 'Linewidth', 2.5);
plot(N_vector, tqp_bvc,'ob', 'MarkerFaceColor', 'b','Linewidth',1.5,'markers',10);

h2 = plot(N_vector, tbuild_bvc, '--b', 'Linewidth', 2.5);
plot(N_vector, tbuild_bvc,'ob', 'MarkerFaceColor', 'b','Linewidth',1.5,'markers',10);

h3 = plot(N_vector, tqp_dmpc, ':r', 'Linewidth', 2.5);
plot(N_vector, tqp_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1.5,'markers',10);

h4 = plot(N_vector, tbuild_dmpc, '--r', 'Linewidth', 2.5);
plot(N_vector, tbuild_dmpc,'or', 'MarkerFaceColor', 'r','Linewidth',1.5,'markers',10);

xlabel('Number of agents');
ylabel('Runtime per Agent [ms]');
[h, icons, plots, s] = legend([h1,h2,h3,h4], 'BVC solve QP', 'BVC build constraint', 'On-demand solve QP', 'On-demand build constraint');
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines([1,5]), 'LineStyle', '-','LineWidth',2); %// modify properties as desired
set(h_lines([3,7]), 'LineStyle', '--','LineWidth',2); %// modify properties as desired
set(gcf,'color','w');










