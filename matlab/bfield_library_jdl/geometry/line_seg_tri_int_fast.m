function [int,p]=line_seg_tri_int_fast(pa,pb,pc,p1,p2)

tol = 1.e-12;
% pa = repmat(pa,nrep,1);
% pb = repmat(pb,nrep,1);
% pc = repmat(pc,nrep,1);
nrep = size(p1,1);
% Calc unit vector normal to plane of Pa-c (gives plane components A-C)
n = cross(pb-pa,pc-pa,2);
n = n./sqrt(sum(n.*n));
n = repmat(n,nrep,1);

% Calculate the position on the line that intersects the plane
mu =  ( dot(n(1,1:3),pa) - dot(n,p1,2) ) ./ dot(n,p2-p1,2);
p = p1 + (p2-p1).*repmat(mu,1,3);

int = zeros(nrep,1);
int(mu >= 0 & mu <= 1,1) = 2;

if any(int == 2)
    %  Determine whether or not the intersection point is bounded by pa,pb,pc
    pa1 = repmat(pa,nrep,1) - p;
    pa1 = pa1./repmat(sqrt(sum(pa1.^2,2)),1,3);
    pa2 = repmat(pb,nrep,1) - p;
    pa2 = pa2./repmat(sqrt(sum(pa2.^2,2)),1,3);
    pa3 = repmat(pc,nrep,1) - p;
    pa3 = pa3./repmat(sqrt(sum(pa3.^2,2)),1,3);
    a1 = dot(pa1,pa2,2);
    a2 = dot(pa2,pa3,2);
    a3 = dot(pa3,pa1,2);
    total = acos(a1) + acos(a2) + acos(a3);
%     abs(total - 2.*pi)
    int(abs(total - 2.*pi) <= tol & int == 2) = 1;
end
int(int ~= 1) =0;

% tol = 1.e-15;
% 
% % Calc unit vector normal to plane of Pa-c (gives plane components A-C)
% n = cross(pb-pa,pc-pa);
% n = n./sqrt(sum(n.*n));
% 
% 
% % Calculate the position on the line that intersects the plane
% mu =  ( dot(n,pa) - dot(n,p1) ) / dot(n,p2-p1);
% p = p1 + mu*(p2-p1);
% 
% if (mu < 0 || mu > 1)    %Intersection not along line segment
%     int = 0;
% else
%     %  Determine whether or not the intersection point is bounded by pa,pb,pc    
%     pa1 = pa - p;
%     pa1 = pa1./sqrt(sum(pa1.^2));
%     pa2 = pb - p;
%     pa2 = pa2./sqrt(sum(pa2.^2));
%     pa3 = pc - p;
%     pa3 = pa3./sqrt(sum(pa3.^2));
%     a1 = dot(pa1,pa2);
%     a2 = dot(pa2,pa3);
%     a3 = dot(pa3,pa1);
%     total = acos(a1) + acos(a2) + acos(a3);
%     if (abs(total - 2.*pi) > tol)
%         int = 0;
%     else
%         int = 1;
%     end
% end