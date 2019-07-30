clc
clear
close all

load('comp_allCA_3.mat')
green = [0.55,0.71,0]; %apple green

% Probability of success plots
prob_bvc = sum(bvc_success, 2) / trials * 100;
prob_bvc2= sum(bvc2_success, 2) / trials * 100;
prob_dmpc = sum(dmpc_success, 2) / trials * 100;
prob_dmpc2 = sum(dmpc2_success, 2) / trials * 100;

figure(1)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([4:4:20]);
ylim([0,105])
xlim([2,22])

h1 = plot(N_vector, prob_bvc,':b', 'Linewidth',2.5);
plot(N_vector, prob_bvc,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);

h2 = plot(N_vector, prob_bvc2,':','Color', green, 'Linewidth',2.5);
plot(N_vector, prob_bvc2,'o','Color', green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);

h3 = plot(N_vector, prob_dmpc,':m', 'Linewidth',2.5);
plot(N_vector, prob_dmpc,'om', 'MarkerFaceColor', 'm','Linewidth',1,'markers',10);

h4 = plot(N_vector, prob_dmpc2,':r','Linewidth',2.5);
plot(N_vector, prob_dmpc2,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);

xlabel('Number of Agents')
ylabel('Success Probability [%]');
[h, icons, plots, s] = legend([h1,h2,h3,h4],{'BVC', 'Soft BVC',...
                             'On-demand (state)', 'On-demand (input)'},'Fontsize',16);
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
set(h,'color','none');


% Completion time
tmean_traj = nanmean(dmpc_trajtime,2);
tstd_traj = nanstd(dmpc_trajtime,1,2);
tmean_traj2 = nanmean(dmpc2_trajtime,2);
tstd_traj2 = nanstd(dmpc2_trajtime,1,2);
tmeanbvc_traj = nanmean(bvc_trajtime,2);
tstdbvc_traj = nanstd(bvc_trajtime,1,2);
tmeanbvc_traj2 = nanmean(bvc2_trajtime,2);
tstdbvc_traj2 = nanstd(bvc2_trajtime,1,2);

figure(2)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([4:4:20]);
xlim([2,22])
ylim([0,25])

h1 = plot(N_vector, tmeanbvc_traj,':b', 'Linewidth',2.5);
plot(N_vector, tmeanbvc_traj,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);

h2 = plot(N_vector, tmeanbvc_traj2,':','Color', green, 'Linewidth',2.5);
plot(N_vector, tmeanbvc_traj2,'o','Color', green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);

h3 = plot(N_vector, tmean_traj,':m', 'Linewidth',2.5);
plot(N_vector, tmean_traj,'om', 'MarkerFaceColor', 'm','Linewidth',1,'markers',10);

h4 = plot(N_vector, tmean_traj2,':r','Linewidth',2.5);
plot(N_vector, tmean_traj2,'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);

xlabel('Number of Agents')
ylabel('Transition Time [s]');
[h, icons, plots, s] = legend([h1,h2,h3,h4],{'BVC', 'Soft BVC',...
                             'On-demand (state)', 'On-demand (input)'},'Fontsize',16);
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
set(h,'color','none');
% set(h, 'Fontsize', 16)


% Plot times
tqp_dmpc = mean(dmpc_tqp, 2);
tqp_dmpc2 = mean(dmpc2_tqp, 2);
tqp_bvc = mean(bvc_tqp, 2);
tqp_bvc2 = mean(bvc2_tqp, 2);


figure(3)
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',20)
hold on;
box on;
xticks([4:4:20]);
xlim([2,22])

h1 = plot(N_vector, tqp_bvc,':b', 'Linewidth',2.5);
plot(N_vector, tqp_bvc,'ob', 'MarkerFaceColor', 'b','Linewidth',1,'markers',10);

h2 = plot(N_vector, tqp_bvc2,':','Color', green, 'Linewidth',2.5);
plot(N_vector, tqp_bvc2,'o','Color', green, 'MarkerFaceColor', green,'Linewidth',1,'markers',10);

h3 = plot(N_vector, tqp_dmpc,':m', 'Linewidth',2.5);
plot(N_vector, tqp_dmpc,'om', 'MarkerFaceColor', 'm','Linewidth',1,'markers',10);

h4 = plot(N_vector, tqp_dmpc2,':r','Linewidth',2.5);
plot(N_vector, tqp_dmpc2, 'or', 'MarkerFaceColor', 'r','Linewidth',1,'markers',10);

xlabel('Number of Agents');
ylabel('Runtime per Agent [ms]');
[h, icons, plots, s] = legend([h1,h2,h3,h4],{'BVC', 'Soft BVC',...
                             'On-demand (state)', 'On-demand (input)'},'Fontsize',16);
h_lines = findobj(icons, 'Type', 'Line');
set(h_lines, 'LineStyle', '-','LineWidth',4); %// modify properties as desired
set(gcf,'color','w');
set(h,'color','none');
