function [pint1,ierr,found_ind1,found_ind2,int_count]=int_curve_curve(line1_r,line1_z,line2_r,line2_z,first,verbose)
% [pint1,ierr,found_ind1,found_ind2,int_count]=int_curve_curve(line1_r,line1_z,line2_r,line2_z,first,verbose)
% Curve is defined by array of points,
% First curve is stepped and both curves linearly interpolated to find intersection with second curve.
% If first == true then first int is returned
% JL 2/2011

if nargin < 5
    first = 0;
end
if nargin < 6
    verbose = 0;
end

int_count = 0;
ierr = 1;
pint1 =  zeros(0,2);
found_ind1 = [];
found_ind2 = [];

for ii1 = 1:length(line1_r) - 1
    p1 = [line1_r(ii1)  ,line1_z(ii1)];
    p2 = [line1_r(ii1+1),line1_z(ii1+1)];

    % Skip zero length
    if p1(1)==p2(1) && p1(2)==p2(2)
        continue
    end

    [pint1_tmp,ierr_tmp,found_ind2_tmp,int_count_tmp]=int_line_curve(p1,p2,line2_r,line2_z,first,verbose);
    
    if ierr_tmp == 0
        idx = int_count + (1:int_count_tmp);
        pint1(idx,1:2)   = pint1_tmp(:,1:2);
        found_ind1(idx)  = ii1;
        found_ind2(idx)  = found_ind2_tmp;        
        int_count = int_count + int_count_tmp;
        ierr = 0;

        if first == 1
            if verbose
                disp('int_curve_curve exiting on first intersection')
            end
            break
        end
    end
end

if int_count == 0
    pint1 = [NaN NaN];
    ierr = 1;
    found_ind1 = NaN;
    found_ind2 = NaN;
end
if int_count > 1 && verbose
    disp('Warning: More than one intersection found')
end

