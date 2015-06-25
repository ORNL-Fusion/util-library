function [pint1,ierr,found_ind1,found_ind2,int_count]=int_curve_curve(line1_r,line1_z,line2_r,line2_z,first,verbose)
% Curve is defined by array of points,
% First curve is stepped and both curves linearly interpolated to find intersection with second curve.
% JL 2/2011

if nargin < 5
    first = 0;
end
if nargin < 6
    verbose = 0;
end

int_count = 0;
found_ind1 = 0;
found_ind2 = 0;
exiting = 0;

for ii1 = 1:length(line1_r) - 1
    p1 = [line1_r(ii1),line1_z(ii1)];
    p2 = [line1_r(ii1+1),line1_z(ii1+1)];
    for ii2 = 1:length(line2_r) - 1
        p3 = [line2_r(ii2),line2_z(ii2)];
        p4 = [line2_r(ii2+1),line2_z(ii2+1)];
        [u1,u2] = int_two_lines(p1,p2,p3,p4);
        dL1 = sqrt(sum((p2-p1).^2));
        dL2 = sqrt(sum((p4-p3).^2));
        
        if ~((sign(dL1) ~= sign(u1)) | (sign(dL2) ~= sign(u2)) | (abs(u1) > 1.d0) | (abs(u2) > 1.d0))
            dR1 = p2(1)-p1(1);
            dZ1 = p2(2)-p1(2);
            theta1 = atan2(dZ1,dR1);
            pint1 = [dL1*u1*cos(theta1)+p1(1),dL1*u1*sin(theta1)+p1(2)];
            
            ierr = 0;
            found_ind1 = ii1;
            found_ind2 = ii2;
            int_count = int_count + 1;
            if first == 1
                if verbose
                    disp('int_line_curve exiting on first intersection')
                end
                exiting = 1;
                break
            end
        end        
    end
    if exiting == 1
        break;
    end
end

if int_count == 0
    pint1 = [0.,0.];
    ierr = 1;
end
if int_count > 1
    disp('Warning: More than one intersection found. Returning last point')
end

