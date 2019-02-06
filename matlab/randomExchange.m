function [po,pf] = randomExchange(N,pmin,pmax,rmin)
max_iter = 200000;
%Generate initial points

pass = false;

while(~pass)
    po = [];
    po(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
    for n = 2:N
        tries = 0;
        pass = false;
        while(~pass && tries <= max_iter)
            candidate = (pmin + (pmax-pmin).*rand(1,3))';
            diff = po - candidate;
            dist = sqrt(sum(diff.^2,1));

            if(dist > rmin)
                po(:,n) = candidate;
                pass = true;   
            end
            tries = tries + 1;
        end
        if (tries > max_iter)
            break;
        end
    end
end

%Make a random permutation of initial states to assign final states
perm = [];
array = 1:N;
array_aux = [];
for i=1:N
    array_aux = [];
    array_aux = array;
    array_aux(array_aux == i) = [];
    
    if i== N
        perm(i) = array;
    elseif (i == N-1 && array_aux(end) == N)
        perm(i) = N;
        array(array == perm(i)) = [];
    else
        j = randi([1 N-i]);
        perm(i) = array_aux(j);
        array(array == array_aux(j)) = [];
    end  
end

for n = 1:N
    pf(:,n) = po(:,perm(n));
end

po = reshape(po,1,3,N);
pf = reshape(pf,1,3,N);