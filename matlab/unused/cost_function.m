% Testing of a strain cost function
clc
close all
clear
err = -1:0.1:1;
vel = -1:0.1:1;
[X,Y] = meshgrid(err,vel);
F = (X.^5)./(-Y+ 0.01);
surf(X,Y,F)
xlabel('Prediction error [m]')
ylabel('Velocity [m/s]')
zlabel('')