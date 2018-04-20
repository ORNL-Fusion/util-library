function [dist,pu,CUTOFF,pu_unchecked,dist_unchecked] = distance_point_to_line_seg(p1,p2,p0)
    % line segment is defined p1 to p2
    % p0 is the reference point
    % endpoints are enforced, see below.
    %
    % pu is point on line seg closest to p0
    % -- If an endpoint is used, CUTOFF = 1
    % in either case pu_unchecked is the closest point on the full line (not limited by end points), similar for dist_unchecked
    
    DEBUG = 0;
    
    nm21 = norm(p2-p1);
    
    dist_unchecked = abs( (p1(2)-p0(2))*(p2(1)-p1(1)) - (p1(1)-p0(1))*(p2(2)-p1(2))  )/nm21;
    u = ( (p0(1)-p1(1))*(p2(1)-p1(1)) + (p0(2)-p1(2))*(p2(2)-p1(2)) ) / nm21^2;
    pu = p1 + u*(p2-p1);
    pu_unchecked = pu;
    
    CUTOFF = 0;
    if u < 0
        pu = p1;
        CUTOFF = 1;
    elseif u > 1
        pu = p2;
        CUTOFF = 1;
    end
    dist = norm(pu-p0);
    
    if DEBUG
        figure; hold on; box on;
        plot([p1(1),p2(1)],[p1(2),p2(2)],'o-')
        plot([p0(1),pu(1)],[p0(2),pu(2)],'kx-')
        plot([p0(1),pu_unchecked(1)],[p0(2),pu_unchecked(2)],'kx-')
    end
end