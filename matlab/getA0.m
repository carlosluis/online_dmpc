function [A0, A0_vel] = getA0(A,K)
Ak = eye(6);
A0 = [];
A0_vel = [];
for  k = 1:K
   Ak = A*Ak;
   A0 = [A0; Ak(1:3,:)];
   A0_vel = [A0_vel; Ak(4:6,:)];
end