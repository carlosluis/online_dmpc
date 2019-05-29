%% %%%%%%%%%%%% 3D VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
[x,y,z] = meshgrid(-5:0.2:5);
figure(1)
colors = distinguishable_colors(N);
% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
while get(gcf,'currentchar')==' '
    for i = 1:N
    h_line(i) = animatedline('LineWidth', 2, 'Color', colors(i,:), 'LineStyle', ':');
    end
    for k = 1:K
        for i = 1:N_cmd
            if k ~= 1
                delete(h_pos(i))
            end
            clearpoints(h_line(i));
            addpoints(h_line(i), hor_rob(1,:,i,k), hor_rob(2,:,i,k), hor_rob(3,:,i,k));
            hold on;
            grid on;
            xlim([phys_limits.pmin(1), phys_limits.pmax(1)])
            ylim([phys_limits.pmin(2), phys_limits.pmax(2)])
            zlim([0, phys_limits.pmax(3)+1.0])
            h_pos(i) = plot3(pos_k_i(1,k,i), pos_k_i(2,k,i), pos_k_i(3,k,i), 'o',...
                             'LineWidth', 2, 'Color',colors(i,:));
            plot3(po(1,1,i), po(1,2,i), po(1,3,i), '^',...
                  'LineWidth', 2, 'Color', colors(i,:));
            plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i), 'x',...
                  'LineWidth', 2, 'Color', colors(i,:));   
        end
        for i = N_cmd + 1: N
            % Plot rouge agents sphere for better visualization
            v = ((x-pos_k_i(1,k,i))/c(i,1)).^(order(i)) + ((y-pos_k_i(2,k,i))/c(i,2)).^(order(i)) + ...
                ((z-pos_k_i(3,k,i))/c(i,3)).^(order(i)) - rmin(i)^order(i);
            patch(isosurface(x,y,z,v,0), 'Facealpha', 0.5, 'FaceColor',...
                           [0.3,0.3,0.3], 'EdgeColor', [0.3,0.3,0.3]);

        end
    drawnow
    end
    pause(0.5)
    clf
end

%% FOR 3D FIGURE
figure(2)
[x,y,z] = meshgrid(-5:0.2:5);
set(gca,'LineWidth',1.2,'TickLength',[0.02 0.02]);
set(gca,'FontSize',14)
for i = 1:N_cmd
    hold on;
    h_pos(i) = plot3(pos_k_i(1,:,i), pos_k_i(2,:,i), pos_k_i(3,:,i),...
                             'LineWidth', 3, 'Color',colors(i,:));
    xlim([phys_limits.pmin(1)-0.2, phys_limits.pmax(1)+0.2])
    ylim([phys_limits.pmin(2)-0.2, phys_limits.pmax(2)+0.2])
    zlim([0, phys_limits.pmax(3)+0.5])
    

end
hold on;
for i = 1:N_cmd
    plot3(po(1,1,i), po(1,2,i), po(1,3,i), '^',...
                  'LineWidth', 5, 'Color', colors(i,:),'Markers',12); 
end

for i = N_cmd + 1: N
            % Plot rouge agents sphere for better visualization
            v = ((x-pos_k_i(1,k,i))/c(i,1)).^(order(i)) + ((y-pos_k_i(2,k,i))/c(i,2)).^(order(i)) + ...
                ((z-pos_k_i(3,k,i))/c(i,3)).^(order(i)) - rmin(i)^order(i);
            patch(isosurface(x,y,z,v,0), 'Facealpha', 0.6, 'FaceColor',...
                           [0.3,0.3,0.3], 'EdgeColor', [0.3,0.3,0.3]);

end
box on;
xticks([-1  1]);
yticks([-1  1]);
zticks([0  2]);
xlabel('x [m]')
ylabel('y [m]');
zlabel('z [m]');

set(gcf,'color','w');

%% FOR 2D FIGURE
figure(2)
[x,y,z] = meshgrid(-1:0.01:1);
set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);
set(gca,'FontSize',20)
hold on;
for i = 1:N_cmd
    h_pos(i) = plot(pos_k_i(1,:,i), pos_k_i(2,:,i),...
                             'LineWidth', 3, 'Color',colors(i,:));
    xlim([-1.2, 1.2])
    ylim([-1.2, 1.2])    

end
for i = 1:N_cmd
    h = plot(po(1,1,i), po(1,2,i), '^',...
                  'LineWidth', 5, 'Color', colors(i,:),'Markers',15); 
end

for i = N_cmd + 1: N
            % Plot rouge agents sphere for better visualization
            v = ((x-pos_k_i(1,k,i))/c(i,1)).^(order(i)) + ((y-pos_k_i(2,k,i))/c(i,2)).^(order(i)) + ...
                ((z-pos_k_i(3,k,i))/c(i,3)).^(order(i)) - rmin(i)^order(i);
            patch(isosurface(x,y,z,v,0), 'Facealpha', 1.0, 'FaceColor',...
                           [0.3,0.3,0.3], 'EdgeColor', [0.3,0.3,0.3]);

end
box on;
% ax = gca;
% ax.LineWidth = 1.2;
xticks([-1  1]);
yticks([-1  1]);
zticks([0  2]);
xlabel('X [m]')
ylabel('Y [m]');
zlabel('Z [m]');

set(gcf,'color','w');
