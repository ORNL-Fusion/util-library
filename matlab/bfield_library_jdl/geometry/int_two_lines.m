function [u1,u2,ierr,pint] = int_two_lines(p1,p2,p3,p4)
% Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
% point is an x-y pair.  Returns two values, which are the normalized distances
% to the intersection point along the line p1 to p2 (u1) and the line
% p3 to p4 (u2).
%
% Equations for line 1 (x1 = [x,y]), and line 2 (x2 = [x,y])
%   _   _        _  _
%   x1 = p1 + u1(p2-p1)
%   _   _        _  _
%   x2 = p3 + u2(p4-p3)
%
% Then solve equastions for u1, u2

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2));
if abs(denom) < eps  % Parallel lines
    u1=NaN;
    u2=NaN;
    pint = NaN;
    ierr = 1;
else     
    u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom;
    if nargout > 1
        u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom;
    end
    if nargout > 3
        pint = p1 + u1*(p2-p1);
    end
    ierr = 0;
end
