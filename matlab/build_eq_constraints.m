function A_eq =  build_eq_constraints(d,l,ndim,deg_poly,T_ctrl_pts)

% From the control points of each derivative of the Bezier curve, we must
% select the first point and the overlapping points between segments
D = zeros(l, 3*(d+1));
for k = 1:l
    if k==1
       D(1,1) = 1; 
    else
       D(k,(k-1)*(d+1)) = 1;
       D(k,(k-1)*(d+1) + 1) = -1;
    end
end

A_eq = augment_array_ndim(D,ndim);

% Make T_ctrl_pts to be a 3D matrix and representing l Bezier segments 
if deg_poly > 0
    T_ctrl_pts_3d{deg_poly} = [];
    D_der{deg_poly} = [];
    for k = 1:deg_poly
        % Select the ctrl points we need and sustract them
        for n = 1:l
            if n == 1
                D_der{k}(1,:) = [T_ctrl_pts{k}(1,:) zeros(1,(l-1)*(d+1))];
            else
                cols = (n-2)*(d+1)+1: (n)*(d+1);
                D_der{k}(n,cols) = [T_ctrl_pts{k}(end,:) -T_ctrl_pts{k}(1,:)];
            end
        end
        % Construct equality constraint matrix by stacking matrices together
        A_eq = [A_eq;augment_array_ndim(D_der{k},3)];
    end
end