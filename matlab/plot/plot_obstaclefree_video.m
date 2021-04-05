%{
Set MATLAB root to 
../data/experiments/20drones_video
before running script
%}
clc
close all
clear all

obstacles_present = 0;

%% Extract data from text files
N = 20;

for i = 1: N
    name = "pred" + num2str(i-1) + ".txt";
    pred(:,:,i) = dlmread(name, '');
end

K = floor(size(pred, 1) / 3);

for i = 1 : N
    l = 1;
    for k = 1 : 3 : 3 * K
        hor(:, :, i, l) = pred(k:k+2, :, i);
        l = l+1;
    end
end

%% %%%%%%%%%%%% 3D VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
colors = distinguishable_colors(N);
[x,y,z] = meshgrid(-5:0.2:5);
pmin = [-2.3,-2.3,0.0];
pmax = [2.3,2.3,2.5];
view(-25,35)
for i = 1:N
h_line(i) = animatedline('LineWidth', 2, 'Color', colors(i,:), 'LineStyle', ':');
end
for k = 1:K
    for i = 1:N
        if k ~= 1
            delete(h_pos(i))
        end
        clearpoints(h_line(i));
        addpoints(h_line(i), hor(1,:,i,k), hor(2,:,i,k), hor(3,:,i,k));
        hold on;
        box on;
        xlim([pmin(1),pmax(1)])
        ylim([pmin(2),pmax(2)])
        zlim([0,pmax(3)])
        h_pos(i) = plot3(hor(1,1,i,k), hor(2,1,i,k), hor(3,1,i,k), 'o',...
                     'LineWidth', 2.0, 'Color',colors(i,:), ...
                     'MarkerEdgeColor','k',...
                     'MarkerFaceColor',colors(i,:),'markers',15);
    end

    if (obstacles_present)
        for i = 1 : N_obs
            % Plot rouge agents sphere for better visualization
            v = ((x - pobs(1,i))/c(i,1)).^(order) + ((y-pobs(2,i))/c(i,2)).^(order) + ...
                ((z-pobs(3,i))/c(i,3)).^(order) - rmin(i)^order;
            patch(isosurface(x,y,z,v,0), 'Facealpha', 0.2, 'FaceColor',...
               [0.3,0.3,0.3], 'EdgeColor', [0.3,0.3,0.3],'LineStyle', 'none');

        end
    end
drawnow
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([-2 0 2])
yticks([-2 0 2])
zticks([0 1 2])
set(gca,'FontSize',16)
set(gcf,'color','w');
F(k) = getframe(gcf) ;
end


%%
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
