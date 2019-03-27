function [po, pf] = random_test_static_rogues(N, N_cmd, pmin, pmax, rmin, E1, order)
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
            diff = E1(:,:,n)*(po - candidate);
            dist = (sum(diff.^order(n),1)).^(1 / order(n));

            if(dist > rmin(n))
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
po = reshape(po, 1, 3, N);

%Generate final points

pass = false;

pf = squeeze(po(:,:,N_cmd + 1:N));
while(~pass)
    if isempty(pf)
        pf(:,1) = (pmin + (pmax-pmin).*rand(1,3))';
    end
    for n = size(pf,2) + 1:N
        tries = 0;
        pass = false;
        while(~pass && tries <= max_iter)
            candidate = (pmin + (pmax-pmin).*rand(1,3))';
            diff = E1(:,:,n)*(pf - candidate);
            dist = (sum(diff.^order(n),1)).^(1/order(n));

            if(dist > rmin(n))
                pf(:,n) = candidate;
                pass = true;   
            end
            tries = tries + 1;
        end
        if (tries > max_iter)
            break;
        end
    end
end
pf = reshape(pf(:,N-N_cmd + 1:N), 1, 3, N_cmd);