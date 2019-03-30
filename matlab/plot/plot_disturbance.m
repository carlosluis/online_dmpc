clc
close all
clear all
load('disturbance.mat')

% Position
figure(1)
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
set(gca,'FontSize',16)
hold on
state = 1;
plot(t, pos_k_i_sample(state,:,1), 'Linewidth', 2.0, 'Color', colors(1,:))
plot(t, ref_sample(state,:,1,1), ':', 'Linewidth', 2.0, 'Color', colors(1,:))
box on;
xlabel("Time [s]");
ylabel("X Position [m]");
ylim([-0.2, 1.4])
yticks([-0.2, 0 , 0.4, 0.8, 1.2])
legend('State', 'Reference');
set(gcf,'color','w');
