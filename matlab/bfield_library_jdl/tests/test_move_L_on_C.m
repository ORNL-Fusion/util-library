clearvars;

% test move_L_on_C with circle
npts = 2001; 
R_C = 0;
Z_C = 0;
radius = 1;
if 1 % nonconstant spacing
    if mod(npts,2)
        theta = linspace(0,pi,(npts+1)/4);
        thetab = linspace(pi,2*pi,npts -(npts+1)/4);        
    else
        theta = linspace(0,pi,npts/4);
        thetab = linspace(pi,2*pi,npts-npts/4);
    end        
    theta = [theta,thetab(2:end)];
else
    theta = linspace(0,2*pi,npts);
end
r = R_C + radius*cos(theta);
z = Z_C + radius*sin(theta);

L = sum(sqrt(diff(r).^2 + diff(z).^2));
circum = 2*pi*radius;
resolve_err = (circum - L)/circum

L_test = L*.73;
[icurve_near_L,err_near_L,R_L,Z_L] = move_L_on_C(L_test,r,z)

dtheta = mod(atan2(Z_L-Z_C,R_L-R_C),2*pi) - theta(1);
dist = dtheta*radius;
err = dist - L_test

figure; hold on; box on;
plot(r,z,'b.-')
plot(R_L,Z_L,'ro')
plot([R_C,r(1)],[Z_C,z(1)],'k-')
plot([R_C,R_L],[Z_C,Z_L],'k-')
title('Distance is CCW')