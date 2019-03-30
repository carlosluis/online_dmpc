clc
clear all;
close all;
load('online_dmpc_1Hz.mat')

%% Get a unified time base for states and commands
% The commands time base is the reference, as it is always
% activated after the full state topic
idx = find((time_fullstate - time_cmds(1,1)) >= 0 & (time_fullstate - time_cmds(1,1)) <= 40);
time_fullstate = time_fullstate(idx);

% Reset time to start at t = 0
time_fullstate = time_fullstate - time_fullstate(1,1);
time_cmds = time_cmds - time_cmds(1,1);
idx_cmds = find(time_cmds <= 40);
time_cmds = time_cmds(idx_cmds);

colors = distinguishable_colors(drone_num);
%%
% Position and commands in X direction
state_to_plot = 1;
figure(1)
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
set(gca,'FontSize',16)
hold on
for i = 1:drone_num
    h_plot(i) =  plot(time_fullstate, fullstate(idx, drones(i), 1, state_to_plot),...
        'LineWidth',2.0,'Color',colors(i,:));
    h_label{i} = ['Drone #' num2str(i)];
    plot(time_cmds(:,i), cmds(idx_cmds, drones(i), state_to_plot),':',...
        'LineWidth',2.0,'Color',colors(i,:));
end
box on;
xlabel("Time [s]");
ylabel("X Position [m]");
legend('State', 'Reference');
set(gcf,'color','w');

% Position and commands in Y direction
state_to_plot = 2;
figure(2)
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
set(gca,'FontSize',16)
hold on
for i = 1:drone_num
    h_plot(i) =  plot(time_fullstate, fullstate(idx, drones(i), 1, state_to_plot),...
        'LineWidth',2.0,'Color',colors(i,:));
    h_label{i} = ['Drone #' num2str(i)];
    plot(time_cmds(:,i), cmds(idx_cmds, drones(i), state_to_plot),':',...
        'LineWidth',2.0,'Color',colors(i,:));
end
box on;
ylim([-1.1 0.6])
xlabel("Time [s]");
ylabel("Y Position [m]");
legend('State', 'Reference');
set(gcf,'color','w');

% Position and commands in Z direction
state_to_plot = 3;
figure(3)
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
set(gca,'FontSize',16)
hold on
for i = 1:drone_num
    h_plot(i) =  plot(time_fullstate, fullstate(idx, drones(i), 1, state_to_plot),...
        'LineWidth',2.0,'Color',colors(i,:));
    h_label{i} = ['Drone #' num2str(i)];
    plot(time_cmds(:,i), cmds(idx_cmds, drones(i), state_to_plot),':',...
        'LineWidth',2.0,'Color',colors(i,:));
end
box on;
xlabel("Time [s]");
ylabel("Z Position [m]");
legend('State', 'Reference');
set(gcf,'color','w');


