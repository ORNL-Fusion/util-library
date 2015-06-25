function surf = plot_vmec_surfaces(wout,phi)
% phi in radians
if nargin < 2
    phi = 0;
end

ntheta = 200;
numsurf = 20;
surf_inds = round(linspace(1,wout.ns,numsurf));

theta = linspace(0,2*pi,ntheta);

for i = 1:length(surf_inds)
    isurf=surf_inds(i);
    for jp = 1:ntheta
        cosmn = cos(wout.xm.*theta(jp)-wout.xn*phi);
        sinmn = sin(wout.xm.*theta(jp)-wout.xn*phi);        
        rsurf(isurf,jp) = sum(wout.rmnc(:,isurf).*cosmn);
        zsurf(isurf,jp) = sum(wout.zmns(:,isurf).*sinmn);
        if wout.asym
            rsurf(isurf,jp) = rsurf(isurf,jp) + sum(wout.rmns(:,isurf).*sinmn);
            zsurf(isurf,jp) = zsurf(isurf,jp) + sum(wout.zmnc(:,isurf).*cosmn);
        end
    end
end

figure; hold on; box on;
plot(rsurf.',zsurf.')

surf.rsurf = rsurf;
surf.zsurf = zsurf;