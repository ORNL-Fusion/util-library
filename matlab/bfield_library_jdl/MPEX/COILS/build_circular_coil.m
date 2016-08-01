function [coil,current] = build_circular_coil(r1,r2,z1,dz,nturns,nlayers,cur,ntheta_per_wind)

if nargin < 8
    ntheta_per_wind = 200;  
end

debug_plots = 0;


% r1 = 1;
% r2 = 1.3;
% z1 = 1;
% dz = 0.3;
% nturns = 3;
% nlayers = 2;
% cur = 1;

% layers in radius, turns in axial
% coils centered around r=0
nwind = nturns*nlayers;
% windings ordered from min Z, min R; R increases fastest
%  r2 -------------------
%     |     |     |     |
%     |  4  |  5  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  2  |  3  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz


% This value includes 0 and 2*pi, so number of current filaments per winding = ntheata_per_wind - 1
% this means number of current filaments per coil = nwind*ntheata_per_wind
% but only nwind*(ntheta_per_wind-1) should carry current
% every ntheta_per_wind filaments should have current = 0





fw = dz/nturns;
fh = (r2-r1)/nlayers;
zwind0 = repmat([z1+fw/2:fw:z1+dz-fw/2],1,nlayers);
rwind0 = repmat([r1+fh/2:fh:r2-fh/2],1,nturns);

thetas = linspace(0,2*pi,ntheta_per_wind);

rwind   = zeros(1,ntheta_per_wind*nwind);
zwind   = zeros(1,ntheta_per_wind*nwind);
current = zeros(ntheta_per_wind*nwind,1);
for iturn = 1:nturns
    for ilayer = 1:nlayers
        iwind = ilayer+(iturn-1)*nlayers;
        ind0 = 1+(iwind-1)*ntheta_per_wind;
        ind1 = ind0 + ntheta_per_wind - 1;
        rwind(ind0:ind1) = rwind0(iwind);
        zwind(ind0:ind1) = zwind0(iwind);
        thetawind(ind0:ind1) = thetas;
        current(ind0:ind1-1) = cur;  % excludes filaments connecting windings
    end
end

coil(:,1) = rwind.*cos(thetawind);
coil(:,2) = rwind.*sin(thetawind);
coil(:,3) = zwind;

if debug_plots
    xwind = rwind.*cos(thetawind);
    ywind = rwind.*sin(thetawind);
    
    figure; hold on; box on;
    plot3(xwind,ywind,zwind,'k-')
    for iwind = 1:nwind
        i0 = 1+(iwind-1)*ntheta_per_wind;
        i1 = i0 + ntheta_per_wind - 1;
        plot3(xwind(i0:i1),ywind(i0:i1),zwind(i0:i1),'b-','linewidth',2)
        if iwind < nwind
            plot3(xwind(i1:i1+1),ywind(i1:i1+1),zwind(i1:i1+1),'r-','linewidth',2)
        end
    end
    
    figure; hold on; box on;
    for i = 1:length(xwind)-1
        if abs(current(i)) < 1e-8
            plot3(xwind(i:i+1),ywind(i:i+1),zwind(i:i+1),'r-')
        else
            plot3(xwind(i:i+1),ywind(i:i+1),zwind(i:i+1),'b-')
        end
    end
end