%{
Data to use with this script:
../data/experiments/10_drone_hulahoop/hoop10_exp.mat
%}
clc
clear all;
close all;
% load ('hoop10_exp.mat');

% Select which mat file you want to load
[FileName,PathName,FilterIndex] = uigetfile('.mat');
load(FileName);

% Transform the status array into something simpler
status = squeeze(status(:,1,:,:));
for i = 1:drone_num
    vehicle_status(:,1,i) = strtrim(join(string(status(:,i,:)),'',2)); 
end

% We want to look only at the CustomControl segment
% We find the beginning and end of 'CustomControl' in vehicle_status
% and find the closest time stamps corresponding to those times

[i,j] = find(contains(vehicle_status(:,1,:),'CustomControl'));
transit_tstart = time_status(min(i));
transit_tstop = time_status(max(i));

% find the start and endpoints of the full_state vector
[val,idx_start] = min(abs(transit_tstart-time_fullstate));
[val,idx_stop] = min(abs(transit_tstop-time_fullstate));

% Reset time to start at t = 0
time_fullstate = time_fullstate(idx_start:idx_stop);

%reset timebase to start at zero
%where zero means the start of the transition
time_fullstate = time_fullstate - time_fullstate(1);

% Make the fullstate a more normal structure, only on the time frame
% required
fullstate = squeeze(fullstate(idx_start:idx_stop,:,1,:));

%% Obtacle definition
N_obs = 4;

% OBSTACLE DEFINITIONS
% first obstacle
rmin(1) = 1.0;
c(1, :) =  [0.35, 1.35, 8.0];

% second obstacle
rmin(2) = 1.0;
c(2, :) =  [0.35, 1.35, 8.0];

% third obstacle
rmin(3) = 1.0;
c(3, :) =  [0.35, 8.0, 1.05];

% fourth obstacle8.1NAN
rmin(4) = 1.0;
c(4, :) =  [0.35, 8.0, 0.45];

% model each obstacle with a seonf order ellipsoid
order = 2;

% Obstacle positions
pobs1 = [0.0, 1.5, 1.0];
pobs2 = [0.0, -1.5, 1.0];
pobs3 = [0.0, 0.0, 0.2];
pobs4 = [0.0, 0.0, 2.0];

pobs = cat(3, pobs1, pobs2, pobs3, pobs4);
pobs = squeeze(pobs);

%% Final positions
pf1 = [1.0, -1.0, 1.0];
pf2 = [1.0, -0.5, 1.0];
pf3 = [1.0, 0.0, 1.0];
pf4 = [1.0, 0.5, 1.0];
pf5 = [1.0, 1.0, 1.0];
pf6 = [-1.0, -1.0, 1.0];
pf7 = [-1.0, -0.5, 1.0];
pf8 = [-1.0, 0.0, 1.0];
pf9 = [-1.0, 0.5, 1.0];
pf10 = [-1.0, 1.0, 1.0];
pf = cat(3, pf1, pf2, pf3, pf4, pf5, pf6, pf7, pf8, pf9, pf10);

%% Dynamic visualization of flight
N = length(drones);
close all
[x,y,z] = meshgrid(-5:0.2:5);
pmin = [-1.5,-1.5,0.0];
pmax = [1.5,1.5,2.0];
tk = 0:0.1:time_fullstate(end);
for i = 1:N
    p(:,:,i) = permute(fullstate(:,i,:), [3,2,1]);
    pk(:,:,i) = spline(time_fullstate,p(:,:,i),tk);
end

% Plotting of drones as boxes
len = 0.15;
width = len;
height = 0.05;
alpha = 0.5;

% plotting hoop as a circle in yz plane\
theta=-pi:0.01:pi;
radius = 0.85 / 2;
center = [0.0, 0.0, 1 + 0.85/2];
y_circle = radius*cos(theta) + center(2);
z_circle = radius*sin(theta) + center(3);
x_circle = center(1)*ones(1, length(y_circle));


X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

colors = distinguishable_colors(N);

% Plot the obstacle and the ellipsoids
figure(1)
for i = 1 : N_obs
    % Plot rouge agents sphere for better visualization
    v = ((x - pobs(1,i))/c(i,1)).^(order) + ((y-pobs(2,i))/c(i,2)).^(order) + ...
        ((z-pobs(3,i))/c(i,3)).^(order) - rmin(i)^order;
    patch(isosurface(x,y,z,v,0), 'Facealpha', 0.3, 'FaceColor',...
                   [0.3,0.3,0.3], 'EdgeColor', [0.3,0.3,0.3],'LineStyle', 'none');
               
    xlim([pmin(1),pmax(1)])
    ylim([pmin(2),pmax(2)])
    zlim([0,pmax(3)])
    hold on
end
box on
plot3(x_circle, y_circle, z_circle, 'k', 'LineWidth', 4)
plot_this = [2, 5, 9, 7];
% plot_this = 1:10;

for i = 1 : N
    plot3(pk(1,1,i), pk(2,1,i),pk(3,1,i), 'o',...
                         'LineWidth', 2.0, 'Color',colors(i,:), ...
                         'MarkerEdgeColor','k',...
                         'MarkerFaceColor',colors(i,:),'markers',20);
                     
    if (ismember(i, plot_this))                 
        plot3(pk(1,:,i), pk(2,:,i), pk(3,:,i),'-', 'Color', colors(i,:), 'Linewidth', 3.0)
    end
end

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([-1 0 1])
yticks([-1 0 1])
zticks([0 1 2])
set(gca,'FontSize',16)
set(gcf,'color','w');


%%
figure(2)
plot3(x_circle, y_circle, z_circle, 'r', 'LineWidth', 3)
hold on
box on
xlim([pmin(1),pmax(1)])
ylim([pmin(2),pmax(2)])
zlim([0,pmax(3)])

plot_this = [2, 5, 9, 8];
% plot_this = 1:10;

for i = 1 : N
    plot3(pk(1,1,i), pk(2,1,i),pk(3,1,i), 'o',...
                         'LineWidth', 0.3, 'Color',colors(i,:), ...
                         'MarkerEdgeColor','k',...
                         'MarkerFaceColor',colors(i,:),'markers',15);
                     
    if (ismember(i, plot_this))                 
        plot3(pk(1,:,i), pk(2,:,i), pk(3,:,i), 'Color', colors(i,:), 'Linewidth', 2)
    end
end

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

%% Plot the min-max distance to goal
l = 1;
dist_goal = [];
for i = 1:drone_num
    diff = p(:,:,i) - repmat(pf(:,:,i), length(time_fullstate), 1)';
    dist_goal(i,:) = sqrt(sum(diff.^2, 1));
end

for i = 1:length(time_fullstate)
        min_dist_goal(i) = min(dist_goal(:,i));
        max_dist_goal(i) = max(dist_goal(:,i));
end

min_dist_goal(min_dist_goal == 0) = NaN;
max_dist_goal(max_dist_goal == 0) = NaN;
% minmin_dgoal = min(min_dist_goal);
% maxmax_dgoal = max(max_dist_goal);

% interpolate for less plotting points
mindgoal_interp = spline(time_fullstate,min_dist_goal(1,:),tk);
maxdgoal_interp = spline(time_fullstate,max_dist_goal(1,:),tk);

% Define colors
light_gray = [0.9, 0.9, 0.9];
light_green = [204,255,204]/255;
figure(3)
set(gca,'LineWidth',1.5,'TickLength',[0.01 0.01]);
set(gca,'FontSize',18)
% Routine for filling areas
success = (0.06)*ones(1,length(tk));
x2 = [tk, fliplr(tk)];
zeros = (0)*ones(1,length(tk));
inBetween = [mindgoal_interp, fliplr(maxdgoal_interp)];
inBetween_success = [zeros, fliplr(success)];
hold on;
box on;
fill(x2, inBetween, light_gray,'LineStyle','none');
fill(x2, inBetween_success,light_green,'LineStyle','none');

h1=plot(tk,maxdgoal_interp,'k', 'Linewidth',1.5);
h2=plot(tk,mindgoal_interp,'k', 'Linewidth',1.5);

plot(tk,success,'--r','LineWidth',1.5);
xlim([0 tk(end)])
xlabel('Time [s]')
ylabel('Distance to Target [m]')
set(gca, 'Layer', 'top');
set(gcf,'color','w');







