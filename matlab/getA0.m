function A0 = getA0(A,K)
Ak = eye(6);
A0 = [];
for  k = 1:K
   Ak = A*Ak;
   A0 = [A0; Ak(1:3,:)];
end