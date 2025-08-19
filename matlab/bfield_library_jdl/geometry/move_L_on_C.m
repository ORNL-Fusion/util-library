function [icurve_near_L,err_near_L,R_L,Z_L] = move_L_on_C(L,rline,zline)
% [icurve_near_L,err_near_L,R_L,Z_L] = move_L_on_C(L,rline,zline)
% Moves a distance L along the curve defined by
% arrays rline,zline.  dL is defined as the
% linear distance between curve points.
% JL 2/2011

if ~isscalar(L)
    error('L must be a scalar')
end
if L < 0
    error('L must be >= 0'); 
end

% Organize inputs
rline = rline(:);
zline = zline(:);
dr = diff(rline);
dz = diff(zline);

segL = hypot(dr,dz);

Ltot = sum(segL);
if L > Ltot
    if L - Ltot <= 10*eps
        R_L = rline(end);
        Z_L = zline(end);
        icurve_near_L = numel(rline);
        err_near_L = 0;
        return;
    else
        error('Requested length %f, exceeds total length %f',L,Ltot)
    end
end

cumSumL = [0;cumsum(segL)];
if L == 0
    ind = 1;
    f = 0;
elseif L == Ltot
    ind = numel(segL);
    f = 1;
else    
    ind = find(cumSumL < L,1,'last');
    if ind == numel(cumSumL)
        ind = numel(segL);
        denom = 0;
    else
        denom = segL(ind);
    end

    if denom == 0
        % Catch zero length segments
        f = 1;
    else
        f = (L - cumSumL(ind))/denom;
    end
end
f = max(0,min(1,f));

R_L = (1-f)*rline(ind) + f*rline(ind+1);
Z_L = (1-f)*zline(ind) + f*zline(ind+1);

% Nearest vertex along curve and error
if f > 0.5
    icurve_near_L = ind + 1;
    err_near_L = (1-f)*segL(ind);
else
    icurve_near_L = ind;
    err_near_L = f*segL(ind);
end
