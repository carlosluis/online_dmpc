function [collisions, neighbours] = check_horizon(hor_k, i, E1, rmin, order)
N = size(hor_k, 2);
collisions = zeros(N, 1);
neighbours = zeros(N, 1);

for j = 1:N %Iterate through the number of obstacles (other agents)
    if(i~=j)
        dist = norm(E1*(hor_k(:,i) - hor_k(:,j)), order);
        collisions(j) = (dist < rmin);
        neighbours(j) = dist < 3*rmin;
    end  
end
