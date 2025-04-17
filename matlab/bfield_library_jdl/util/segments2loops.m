function [r2,z2] = segments2loops(r,z,tol)
% r,z : 2‑by‑N arrays of segment endpoints
% r2,z2 : 1‑by‑m cell arrays of ordered loop coordinates
if nargin<3, tol = 1e-12; end
N = size(r,2); left = true(1,N);
r2 = {}; z2 = {};
while any(left)                              % build each loop
    k = find(left,1);                        % start with a free segment
    R = r(:,k);  Z = z(:,k);  left(k) = false;
    grown = true;
    while grown                              % extend head/tail
        grown = false;
        for j = find(left)                   % test unused segment j
            if  norm([r(1,j);z(1,j)]-[R(end);Z(end)])<tol
                R = [R; r(2,j)]; Z = [Z; z(2,j)];
            elseif norm([r(2,j);z(2,j)]-[R(end);Z(end)])<tol
                R = [R; r(1,j)]; Z = [Z; z(1,j)];
            elseif norm([r(1,j);z(1,j)]-[R(1);Z(1)])<tol
                R = [r(2,j); R]; Z = [z(2,j); Z];
            elseif norm([r(2,j);z(2,j)]-[R(1);Z(1)])<tol
                R = [r(1,j); R]; Z = [z(1,j); Z];
            else
                continue
            end
            left(j) = false;                 % segment consumed
            grown = true;
            break
        end
    end
    if norm([R(1);Z(1)]-[R(end);Z(end)])>tol % close if not closed
        R(end+1) = R(1); Z(end+1) = Z(1);
    end
    r2{end+1} = R;                           % store loop
    z2{end+1} = Z;
end
end
