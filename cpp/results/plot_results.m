clc
close all
clear
view_states = false;
view_animation = true;

M = dlmread('trajectories.txt','');
N = M(1,1);
N_cmd = M(1,2);
K = size(M,2);
pmin = M(1,3:5);
pmax = M(1,6:8);

po = M(2:4,1:N);
po = reshape(po,1,3,N);
pf = M(5:7,1:N_cmd);
pf = reshape(pf,1,3,N_cmd);
%%
start = 8;
final = start + 3*N_cmd-1;
all_pos = M(start:final,:);
pk = [];

for i=1:N_cmd
    pk(:,:,i) = all_pos(3*(i-1)+1:3*i,:);
end

T = 0.01*(size(pk,2)-1);
t = 0:0.01:T; 

if view_states
    for i = 1:N_cmd  
        figure(1)
        diff = pk(:,:,i) - repmat(pf(:,:,i),length(t),1)';
        dist = sqrt(sum(diff.^2,1));
        plot(t, dist, 'LineWidth',1.5);
        grid on;
        hold on;
        xlabel('t [s]')
        ylabel('Distance to target [m]');

        figure(2)
        subplot(3,1,1)
        plot(t,pk(1,:,i),'LineWidth',1.5);
        hold on;
        plot(t,pmin(1)*ones(length(t),1),'--r','LineWidth',1.5);
        plot(t,pmax(1)*ones(length(t),1),'--r','LineWidth',1.5);
        ylabel('x [m]')
        xlabel ('t [s]')
        grid on;

        subplot(3,1,2)
        plot(t,pk(2,:,i),'LineWidth',1.5);
        hold on;
        plot(t,pmin(2)*ones(length(t),1),'--r','LineWidth',1.5);
        plot(t,pmax(2)*ones(length(t),1),'--r','LineWidth',1.5);
        ylabel('y [m]')
        xlabel ('t [s]')
        grid on;


        subplot(3,1,3)
        plot(t,pk(3,:,i),'LineWidth',1.5);
        hold on;
        plot(t,pmin(3)*ones(length(t),1),'--r','LineWidth',1.5);
        plot(t,pmax(3)*ones(length(t),1),'--r','LineWidth',1.5);
        ylabel('z [m]')
        xlabel ('t [s]')
        grid on;

    end
end

% Downsample for better visualization
pk = pk(:,1:20:end,:);

%% Animation of transition
if view_animation
    figure(3)
    colors = distinguishable_colors(N);
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'currentchar',' ')
    while get(gcf,'currentchar')==' '
        for k = 1:size(pk, 2)
            for i = 1:N
                hold on;
                grid on;
                xlim([pmin(1),pmax(1)])
                ylim([pmin(2),pmax(2)])
                zlim([0,pmax(3)])
                if i <= N_cmd

                    plot3(pk(1,k,i),pk(2,k,i),pk(3,k,i),'o',...
                        'LineWidth',2,'Color',colors(i,:));
                    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                          'LineWidth',2,'Color',colors(i,:));
                    plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'x',...
                          'LineWidth',2,'Color',colors(i,:)); 
                else
                    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'^',...
                          'LineWidth',2,'Color',colors(i,:));
                end

            end
        drawnow
        end
        clf
        pause(0.1)
    end
end


