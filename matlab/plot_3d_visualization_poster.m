%% %%%%%%%%%%%% 3D VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
[x,y,z] = meshgrid(-5:0.1:5);
figure(1)
colors = distinguishable_colors(N);
colors(1,:) = [0,0.73,0.98];
colors(3,:) = [0.55,0.75,0.11];
colors(4,:) = [0.66,0.66,0.66];
% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'currentchar',' ')
set(0,'DefaultFigureColor',[1 1 1])
% set(gcf,'Renderer','opengl')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
set(gca,'FontSize',20)
    for i = 1:N
    h_line(i) = animatedline('LineWidth',4,'Color',colors(i,:),'LineStyle',':','markers',3);
    end
for k = 1:K
    for i = 1:N_cmd
        if k ~= 1
%                 delete(h_pos(i))
        end
        clearpoints(h_line(i));
        addpoints(h_line(i), hor_rob(1,:,i,k), hor_rob(2,:,i,k), hor_rob(3,:,i,k));
        hold on;
        box on;
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('z [m]')
        ax = gca;
        ax.LineWidth = 5;
        xticks([-1  1]);
        yticks([-1  1]);
        xlim([phys_limits.pmin(1), phys_limits.pmax(1)])
        ylim([phys_limits.pmin(2), phys_limits.pmax(2)])
        zlim([-10, 10])
        plot(pos_k_i(1,k,i), pos_k_i(2,k,i), 'o',...
                         'LineWidth', 1.5, 'Color',colors(i,:), ...
                         'MarkerEdgeColor','k',...
                         'MarkerFaceColor',colors(i,:),'markers',25);  
    end
    if k==1
        xh = get(gca,'xlabel'); % handle to the label object
        p = get(xh,'position'); % get the current position property
        p(2) = p(2)/1.2 ;        % double the distance, 
                               % negative values put the label below the axis
        set(xh,'position',p)   % set the new position
        yh = get(gca,'ylabel'); % handle to the label object
        p = get(yh,'position'); % get the current position property
        p(1) = p(1)/1.1 ;        % double the distance, 
                               % negative values put the label below the axis
        set(yh,'position',p)   % set the new position
    end
    for i = N_cmd + 1: N + N_obs
        % Plot rouge agents sphere for better visualization
        v = ((x-pos_k_i(1,k,i))/c(i,1)).^(order(i)) + ((y-pos_k_i(2,k,i))/c(i,2)).^(order(i)) + ...
            ((z-pos_k_i(3,k,i))/c(i,3)).^(order(i)) - (rmin(i)-0.2)^order(i);
        patch(isosurface(x,y,z,v,0), 'Facealpha', 1.0, 'FaceColor',...
                       [0.1,0.1,0.1], 'EdgeColor', [0.1,0.1,0.1]);

    end
    for i=1:N_cmd
    plot3(po(1,1,i), po(1,2,i), po(1,3,i),'o',...
                  'LineWidth',2,'Color',colors(i,:),...
                  'MarkerEdgeColor','k',...
                  'MarkerFaceColor',colors(i,:),'markers',25);
    end
    for i =1:N_cmd
    plot3(pf(1,1,i), pf(1,2,i), pf(1,3,i),'d',...
                  'LineWidth',2,'Color',colors(i,:),...
                  'MarkerEdgeColor','k',...
                  'MarkerFaceColor',colors(i,:),'markers',15);
    end
    F(k) = getframe(gcf) ;
    drawnow
%         pause(h)
    set(gca,'color','none');
end

% create the video writer with 1 fps
writerObj = VideoWriter('myVideo.avi', 'Uncompressed AVI');
writerObj.FrameRate = 5;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


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
xlabel('X [m]')
ylabel('Y [m]');
zlabel('Z [m]');

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
