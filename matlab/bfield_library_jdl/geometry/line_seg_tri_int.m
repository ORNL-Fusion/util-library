function [int,p]=line_seg_tri_int(pa,pb,pc,p1,p2)

tol = 1.e-12;

% Calc unit vector normal to plane of Pa-c (gives plane components A-C)
n(1) = (pb(2) - pa(2))*(pc(3) - pa(3)) - (pb(3) - pa(3))*(pc(2) - pa(2));
n(2) = (pb(3) - pa(3))*(pc(1) - pa(1)) - (pb(1) - pa(1))*(pc(3) - pa(3));
n(3) = (pb(1) - pa(1))*(pc(2) - pa(2)) - (pb(2) - pa(2))*(pc(1) - pa(1));
Ln=sqrt(sum(n.^2));
n=n./Ln;

% Calculate plane component D
d = - n(1)*pa(1) - n(2)*pa(2) - n(3)*pa(3);

% Calculate the position on the line that intersects the plane
denom = n(1)*(p2(1) - p1(1)) + n(2)*(p2(2) - p1(2)) + n(3)*(p2(3) - p1(3));

if (abs(denom) < tol)
    int = 0;
    error('should never get here')
else
    
    mu = - (d + n(1) * p1(1) + n(2) * p1(2) + n(3) * p1(3)) / denom;
    p(1) = p1(1) + mu * (p2(1) - p1(1));
    p(2) = p1(2) + mu * (p2(2) - p1(2));
    p(3) = p1(3) + mu * (p2(3) - p1(3));
        
    if (mu < 0 || mu > 1)    %Intersection not along line segment
        int = 0;
    else        
        %  Determine whether or not the intersection point is bounded by pa,pb,pc
        pa1(1) = pa(1) - p(1);
        pa1(2) = pa(2) - p(2);
        pa1(3) = pa(3) - p(3);
        pa1 = pa1./sqrt(sum(pa1.^2));
        pa2(1) = pb(1) - p(1);
        pa2(2) = pb(2) - p(2);
        pa2(3) = pb(3) - p(3);
        pa2 = pa2./sqrt(sum(pa2.^2));
        pa3(1) = pc(1) - p(1);
        pa3(2) = pc(2) - p(2);
        pa3(3) = pc(3) - p(3);
        pa3 = pa3./sqrt(sum(pa3.^2));
        a1 = pa1(1)*pa2(1) + pa1(2)*pa2(2) + pa1(3)*pa2(3);
        a2 = pa2(1)*pa3(1) + pa2(2)*pa3(2) + pa2(3)*pa3(3);
        a3 = pa3(1)*pa1(1) + pa3(2)*pa1(2) + pa3(3)*pa1(3);
        total = (acos(a1) + acos(a2) + acos(a3));
%         abs(total - 2.*pi)
        if (abs(total - 2.*pi) > tol)
            int = 0;
        else
            int=1;
        end
    end
end

end