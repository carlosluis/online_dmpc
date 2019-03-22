clc
close all
clear
% simulation of a rouge agent moving at constant speed

speed = 0.5; % constant speed of the agent in m/s
T = 20;
Ts = 0.2;
t = 0:Ts:T;

max_angle = 30*pi/180; % degrees

p = zeros(3, length(t));
v = zeros(3, length(t));

p(:,1) = [0.0 ; 0.0; 1.0]; % initial position

% generate a random initial orientation for the velocity vector
euler = rnd_yrp(pi);

v(:,1) = eul2rotm(euler') * [0.0; 0.0; speed];
count  = 1;
for k = 2:length(t)
    % generate random rotation
    v(:,k) =  eul2rotm(euler') * v(:,k - 1);
    p(:,k) = p(:, k -1) + Ts * v(:,k);
end

N = 1;
figure(1)
colors = distinguishable_colors(N);
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
    for k = 1:length(t)
        for i = 1:N
            xlim([-2.0,2.0])
            ylim([-2.0, 2.0])
            zlim([0, 2.0])
            h_pos(i) = plot3(p(1,k,i), p(2,k,i), p(3,k,i), 'o',...
                             'LineWidth', 2, 'Color',colors(i,:));
            hold on;
            grid on;

        end
    drawnow
    end
    pause(0.5)
    clf
end
