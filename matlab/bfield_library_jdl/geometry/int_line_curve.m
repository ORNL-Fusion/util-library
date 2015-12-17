function [pint1,ierr,found_ind,int_count]=int_line_curve(p1,p2,line_r,line_z,first,verbose)
% Curve is defined by array of points,
% linearly interpolated to find intersection with line.
% JL 2/2011

if nargin < 5
    first = 0;
end
if nargin < 6
    verbose = 0;
end

int_count(1) = 0;
found_ind(1) = 0;


nn = length(line_r)-1;
for ii = 1:nn
    p3 = [line_r(ii),line_z(ii)];
    p4 = [line_r(ii+1),line_z(ii+1)];
    [u1,u2,ierr2] = int_two_lines(p1,p2,p3,p4);

    if ierr2 ~= 0  % Parallel line segments
        continue; 
    end
    
    if ii == nn 
        test = u1 >= 0 && u1 <= 1 && u2 >= 0 && u2 <= 1;
    else
        test = u1 >= 0 && u1 <= 1 && u2 >= 0 && u2 < 1;
    end
    
    if test
        int_count = int_count + 1;
        pint1(int_count,:) = p1 + u1*(p2-p1);        
        ierr = 0;
        found_ind(int_count) = ii;        
        if first == 1
            if verbose
                disp('int_line_curve exiting on first intersection')
            end
            break
        end
    end
end

if int_count > 1 && verbose
    disp('Warning: More than one intersection found.')
end

if int_count(1) == 0
    pint1(1,1:2) = NaN;    
    found_ind(1) = NaN;
    ierr = 1;
end





