function [Rax,Zax,dist] = find_axis(bfield,phistart,R0,Z0)
% bfield should contain nsym
tic;
if nargin < 3
    rcoil = sqrt(bfield.coil(:,1).^2 + bfield.coil(:,2).^2);
    zcoil = bfield.coil(:,3);
    R0 = mean(rcoil);
    Z0 = mean(zcoil);
end

X0=[R0,Z0];
OPTIONS=optimset('tolfun',1e-9,'tolx',1e-10,'Display','off');
XLB=[];
XUB=[];


% global params for getdistsq
nowarn = 1;
dphi = 0.5*pi/180;
rnst = 2*pi/bfield.nsym/dphi;
nsteps = round(rnst);
if abs(nsteps - rnst) > 1e-6
    error('Make this divisible')
end


[X]=lsqnonlin(@getdistsq,X0,XLB,XUB,OPTIONS);

[dist_sq,Rax,Zax]=getdistsq(X);
dist=sqrt(dist_sq(1)^2 + dist_sq(2)^2);
fprintf('Found axis at [R,Z] = [%f, %f], err = %e \n',Rax,Zax,dist);
fprintf('Axis search took %f seconds.\n',toc); tic;


    function [dist_sq,r,z]=getdistsq(X)
        
        
        f = follow_fieldlines_rzphi_dphi(bfield,X(1),X(2),phistart,dphi,nsteps,nowarn);
        
        dist_sq=[f.r(1)-f.r(end),f.z(1)-f.z(end)];
        if nargout > 1
            r = f.r(end);
            z = f.z(end);
        end
    end
end


