function [rr,dd] = plot_IR_data_raw(plotit);
% clearvars;
if nargin < 1
    plotit = 1;
end

d = load('C:\Work\MPEX\shot 7445.mat');
% d = load('C:\Work\MPEX\shot 6547.mat');

if 0
    figure; hold on; box on;
    h= pcolor(d.Frame);
    hc = colorbar;
    set(h,'edgecolor','none')
end

nw = size(d.Frame,2);
nh = size(d.Frame,1);

px_per_cm = 12.146;
dr = nw/px_per_cm;
dz = nh/px_per_cm;

r = linspace(0,dr,nw+1);
z = linspace(0,dz,nh+1);
data = flipud(d.Frame);


z0 = 17.2;
r0 = 25.75;

icount = 1;
for i = 1:length(r)-1
    for j = 1:length(z)-1
        rcell(:,icount) = [r(i),r(i+1),r(i+1),r(i)];
        zcell(:,icount) = [z(j),z(j),z(j+1),z(j+1)];
        dcell(:,icount) = data(j,i);
        icount = icount + 1;
    end
end
% r = linspace(0,nw/px_per_cm);
% z = linspace(0,nh)/px_per_cm;

if plotit
    figure;hold on; box on;
    patch(rcell-r0,zcell-z0,dcell,'edgecolor','none')
    plot([0,0],[-z0,z(end)-z0],'k-')
    plot([-r0,r(end)-r0],[0,0],'k-')
    % h= pcolor(r,z,flipud(d.Frame));
    hc = colorbar;
    % set(h,'edgecolor','none')
end

theta = linspace(0,2*pi,100);

r1 = 2.5;
xx = r1*cos(theta);
yy = r1*sin(theta);
if plotit
    plot(xx,yy,'k')
end

r1 = 1.8;
xx = r1*cos(theta);
yy = r1*sin(theta);
if plotit
    plot(xx,yy,'k')
end

theta1 = 0.7*2*pi;
ninterp = 100;
rr = linspace(0,5,ninterp);
xx = rr.*cos(theta1);
yy = rr.*sin(theta1);
if plotit
    plot(xx,yy,'k')
    xlabel('X [cm]')
    ylabel('\DeltaT')
end

% Interpolate along line
rinterp = linspace(0,dr,nw)-r0;
zinterp = linspace(0,dz,nh)-z0;
dd = interp2(rinterp,zinterp,data,xx,yy);
if plotit
    figure; hold on; box on
    plot(rr,dd)
    xlabel('r [cm]')
    ylabel('\DeltaT')
end

