function [A,B] = getABmodel(h)

% Define model parameters for the quad + controller system
zeta_xy = 0.6502;
tau_xy = 0.3815;
omega_xy = 1/tau_xy;
zeta_z = 0.9103;
tau_z = 0.3;
omega_z = 1/tau_z;

% Matrices of the discrete state space model
% State = [x y z vx vy vz]
A = [1  0 0 h 0 0;
     0 1 0 0 h 0;
     0 0 1 0 0 h;
     -h*omega_xy^2 0 0 -2*omega_xy*h*zeta_xy 0 0;
     0 -h*omega_xy^2 0 0 -2*omega_xy*h*zeta_xy 0;
     0 0 -h*omega_z^2 0 0 -2*omega_z*h*zeta_z];
 
B = [zeros(3,3);
     h*omega_xy^2 0 0;
     0 h*omega_xy^2 0;
     0 0 h*omega_z^2]; 