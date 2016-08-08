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
if L < eps
    ind = 1;
else
    ind = find(SumL < L,1,'last');
end
f = (L - SumL(ind))/dL(ind+1);
R_L = f*(rline(ind+1)-rline(ind)) + rline(ind);
Z_L = f*(zline(ind+1)-zline(ind)) + zline(ind);
if f > 0.5
    icurve_near_L = ind+1;
    err_near_L=(1-f)*dL(ind+1);
else
    icurve_near_L = ind;
    err_near_L=f*dL(ind+1);
end
