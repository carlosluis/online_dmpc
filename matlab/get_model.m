function [model, inv_model] = get_model(h, model_params)

% This is bad practice, but here we're doing it a safe way for a specific purpose
% Unpack the structure model_params

names = fieldnames(model_params);
for i = 1:length(names)
   eval([names{i} '=model_params.' names{i} ';']); 
end

% Matrices of the discrete state space model
% State = [x y z vx vy vz]
model.A = [1  0 0 h 0 0;
           0 1 0 0 h 0;
           0 0 1 0 0 h;
           -h*omega_xy^2 0 0 1-2*omega_xy*h*zeta_xy 0 0;
           0 -h*omega_xy^2 0 0 1-2*omega_xy*h*zeta_xy 0;
           0 0 -h*omega_z^2 0 0 1-2*omega_z*h*zeta_z];
 
model.B = [zeros(3, 3);
           h*omega_xy^2 0 0;
           0 h*omega_xy^2 0;
           0 0 h*omega_z^2]; 
 
inv_model.A = [1 0 0 -(1-2*omega_xy*h*zeta_xy)/(h*omega_xy^2) 0 0;
               0 1 0 0 -(1-2*omega_xy*h*zeta_xy)/(h*omega_xy^2) 0;
               0 0 1 0 0 -(1-2*omega_z*h*zeta_z)/(h*omega_z^2)];
  
inv_model.B = diag([1/(h*omega_xy^2), 1/(h*omega_xy^2), 1/(h*omega_z^2)]);