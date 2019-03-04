% Testing of a strain cost function
clc
close all
clear
err = -2:0.05:2;
vel = -2:0.05:2;
[X,Y] = meshgrid(err,vel);
F = (X.^2)./(Y+0.01);
surf(X,Y,F)
xlabel('Error [m]')
ylabel('Velocity [m/s]')
zlabel('Cost function')