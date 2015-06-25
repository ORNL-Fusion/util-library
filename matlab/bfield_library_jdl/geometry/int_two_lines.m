function [u1,u2] = int_two_lines(p1,p2,p3,p4)
% ; Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
% ; point is an x-y pair.  Returns two values, which are the distance
% ; between p1 and p2 of the intersection, and between p3 and p4.  This
% ; distance is normalized to the length of the lines.

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2));
if denom == 0 
    u1=1e30;
    u2=1e30;
else 
    u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom;
    u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom;
end
