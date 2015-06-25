function [icurve_near_L,err_near_L,R_L,Z_L] = move_L_on_C(L,rline,zline)
% Moves a distance L along the curve defined by
% arrays rline,zline.  dL is defined as the
% linear distance between curve points.
% JL 2/2011

nC = length(rline);
dL = zeros(1,nC);
dL(2:nC) = sqrt( (rline(1:nC-1)-rline(2:nC)).^2 + (zline(1:nC-1)-zline(2:nC)).^2 );

Ltot = sum(dL);
if Ltot < L
    fprintf('Requested length %f, total length %f\n',[L,Ltot])
    error('Error: Ltot < L')
end

SumL = cumsum(dL);
delt = L-SumL;
[~,ind_diff] = min(abs(delt));
diff = delt(ind_diff);

if diff > 0
    dL_step = dL(ind_diff+1);
    ii = ind_diff;
    theta = atan2(zline(ii+1)-zline(ii),rline(ii+1)-rline(ii));
    if dL_step == 0
        error(' here')
    end
    R_L = rline(ii) + diff*cos(theta);
    Z_L = zline(ii) + diff*sin(theta);
elseif diff < 0
    dL_step = dL(ind_diff);
    ii = ind_diff;
    theta = atan2(zline(ii-1)-zline(ii),rline(ii-1)-rline(ii));
    if dL_step == 0
        error(' here')
    end
    R_L = rline(ii) - diff*cos(theta);
    Z_L = zline(ii) - diff*sin(theta);
else
    R_L = rline(ind_diff);
    Z_L = zline(ind_diff);
end

icurve_near_L = ind_diff;
err_near_L = diff;

