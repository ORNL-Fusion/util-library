function [pint1,ierr,found_ind,int_count]=int_line_curve(p1,p2,line_r,line_z,first,verbose)
% ; Curve is defined by array of points,
% ; linearly interpolated to find intersection with line.
% ; JL 2/2011

if nargin < 5
    first = 0;
end
if nargin < 6
    verbose = 0;
end

int_count = 0;
found_ind = 0;

for ii = 1:length(line_r)-1
    p3 = [line_r(ii),line_z(ii)];
    p4 = [line_r(ii+1),line_z(ii+1)];
    [u1,u2] = int_two_lines(p1,p2,p3,p4);
    dL1 = sqrt(sum((p2-p1).^2));
    dL2 = sqrt(sum((p4-p3).^2));
        
        if (sign(dL1) ~= sign(u1)) | (sign(dL2) ~= sign(u2)) | (abs(u1) > 1.d0) | (abs(u2) > 1.d0)
        else
            dR1 = p2(1)-p1(1);
            dZ1 = p2(2)-p1(2);
            theta1 = atan2(dZ1,dR1);
            pint1 = [dL1*u1*cos(theta1)+p1(1),dL1*u1*sin(theta1)+p1(2)];
            
%             dR2 = p4(1)-p3(1);
%             dZ2 = p4(2)-p3(2);
%             theta2 = atan2(dZ2,dR2);
%             pint2 = [dL2*u2*cos(theta2)+p3(1),dL2*u2*sin(theta2)+p3(2)];
            ierr = 0;
            found_ind = ii;
            int_count = int_count + 1;
            if first == 1
                if verbose
                    disp('int_line_curve exiting on first intersection')
                end
                break
            end
        end

end

if int_count == 0
    pint1 = [0.,0.];
%     pint2 = 0.;
    ierr = 1;
end
if int_count > 1 & verbose
    disp('Warning: More than one intersection found. Returning last point')
end





