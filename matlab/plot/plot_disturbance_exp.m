%{
Data to be used in this script:
- ../data/experiments/continuous_replanning_1Hz.mat
- ../data/experiments/event_replanning.mat
%}
clc
clear all;
close all;
load ('event_replanning.mat');

% Select which mat file you want to load
% [FileName,PathName,FilterIndex] = uigetfile('.mat');
% load(FileName);

%% Get a unified time base for states and commands
% The commands time base is the reference, as it is always
% activated after the full state topic
idx = find((time_fullstate - time_cmds(1,1)) >= 0);
time_fullstate = time_fullstate(idx);

% Reset time to start at t = 0
time_fullstate = time_fullstate - time_fullstate(1,1);
time_cmds = time_cmds - time_cmds(1,1);

colors = distinguishable_colors(drone_num);
drone_num = 1;
%%
% Position and commands in X direction
state_to_plot = 2;
figure(1)
for i = 1:drone_num
    h_plot(i) =  plot(time_fullstate, fullstate(idx, drones(i), 1, state_to_plot),...
        'LineWidth',2.0,'Color',colors(i,:));
    hold on;
    plot(time_cmds(:,i), cmds(:, drones(i), state_to_plot),':',...
        'LineWidth',1.5,'Color',colors(i,:));
    
end
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
set(gca,'FontSize',16)
xlabel("Time [s]");
xlim([0 84])
ylim([-2, 2])
ylabel("Y Position [m]");
legend('State', 'Reference')
set(gcf,'color','w');