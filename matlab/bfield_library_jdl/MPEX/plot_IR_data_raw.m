% function [rr,dd] = plot_IR_data_raw(shot,plotit)
function plot_IR_data_raw
% if nargin < 2
%     plotit = 1;
% end
clearvars;
% shot = 7445;
shot = 7410;
plotit = 1;

data_path = 'C:\Work\MPEX\Camera\';

files = dir(data_path);

icount = 0;
for i = 3:length(files)
    itmp = strfind(files(i).name,num2str(shot));
    if ~isempty(itmp)
        ishot = i;
        icount = icount + 1;   
    end
end
if icount == 0
    error(['Could not find camera data for shot ',num2str(shot)])
end
if icount > 1
    error(['Found multiple files matching shot: ',num2str(shot)])
end

fprintf('Using file %s \n',files(ishot).name)

d = load(strcat(data_path,files(ishot).name));
if 0  % plot raw data
    figure; hold on; box on;
    h= pcolor(d.Frame);
    hc = colorbar;
    set(h,'edgecolor','none')
end

nw = size(d.Frame,2);
nh = size(d.Frame,1);

px_per_cm = 12.146;
dx = nw/px_per_cm;
dy = nh/px_per_cm;

x = linspace(0,dx,nw+1);
y = linspace(0,dy,nh+1);
data = flipud(d.Frame);

% % 
% y0 = 17.2;
% x0 = 25.75;

% y0 = 18;
% x0 = 26.5;

x0_guess = x(end)/2;
y0_guess = y(end)/2-2.5;
% z0_guess = y(end)/2;


x0 = x0_guess;
y0 = y0_guess;

icount = 1;
xcell = zeros(4,nh*nw);
ycell = zeros(4,nh*nw);
dcell = zeros(1,nh*nw);
for i = 1:nw
    for j = 1:nh
        xcell(:,icount) = [x(i),x(i+1),x(i+1),x(i)];
        ycell(:,icount) = [y(j),y(j),y(j+1),y(j+1)];
        dcell(icount) = data(j,i);
        icount = icount + 1;
    end
end
xmean = mean(xcell);
ymean = mean(ycell);


if plotit
    figure;hold on; box on;
    patch(xcell-x0,ycell-y0,dcell,'edgecolor','none')
    plot([0,0],[-y0,y(end)-y0],'k-')
    plot([-x0,x(end)-x0],[0,0],'k-')
    colorbar;
    xlabel('X [cm]','fontsize',14)
    title('\DeltaT','fontsize',14)
    ylabel('Y [cm]','fontsize',14)

    % x = linspace(0,nw/px_per_cm);
    % y = linspace(0,nh)/px_per_cm;
    % h= pcolor(x,y,flipud(d.Frame));
    % set(h,'edgecolor','none')    
end



theta = linspace(0,2*pi,100);

verbose = 1;
tolx = 1e-6;
tolfun = 1e-6;

if verbose == 0
    qval = 'off';
elseif verbose == 1
    qval = 'final';
elseif verbose == 2
    qval = 'iter';
else
    qval = 'final';
end
opts=optimoptions('lsqnonlin','TolFun',tolfun,'TolX',tolx,'Display',qval,'Algorithm','levenberg-marquardt','TypicalX',[100000000,100000000]);
opts.TolFun=tolfun;
opts.TolX=tolx;
opts.Algorithm = 'levenberg-marquardt';

radius_eval = 2.;  %cm
revals = linspace(0,radius_eval,30);
tevals = linspace(0,2*pi,40); tevals(end) = [];
[rr,tt] = meshgrid(revals,tevals);
xevals = rr.*cos(tt);
yevals = rr.*sin(tt);
xinterp = linspace(0,dx,nw)-x0;
yinterp = linspace(0,dy,nh)-y0;

x00=[1,0];
% xub = [max(x)-x0,max(y)-y0]-RTARG;
% xlb = [min(x)-x0,min(y)-y0]+RTARG;
% xfinal=lsqnonlin(@minfun,x00,xlb,xub,opts);
xfinal=lsqnonlin(@minfun,x00,[],[],opts);

% plot(XTARG+xfinal(1),YTARG+xfinal(2),'m')
% plot(XTARG+x00(1),YTARG+x00(2),'y--')
plot(xfinal(1),xfinal(2),'mx')
plot(x00(1),x00(2),'mo')


% plot(
xfinal

% return;

x0 = x0 + xfinal(1);
y0 = y0 + xfinal(2);

figure;hold on; box on;
% clf; hold on; box on;
patch(xcell-x0,ycell-y0,dcell,'edgecolor','none')
plot([0,0],[-y0,y(end)-y0],'k-')
plot([-x0,x(end)-x0],[0,0],'k-')
colorbar;
xlabel('X [cm]','fontsize',14)
title('\DeltaT','fontsize',14)
ylabel('Y [cm]','fontsize',14)



r1 = radius_eval;
xx = r1*cos(theta);
yy = r1*sin(theta);
if plotit
    plot(xx,yy,'m')
end

r1 = 1.8;
xx = r1*cos(theta);
yy = r1*sin(theta);
if plotit
    plot(xx,yy,'k')
end

xeval = radius_eval*cos(theta);
yeval = radius_eval*sin(theta);
isin = inpolygon(xmean-x0,ymean-y0,xeval,yeval);
[dmax,imax] = max(dcell(isin));
xtmp =  xmean(isin) - x0;
ytmp =  ymean(isin) - y0;

theta1 = atan2(ytmp(imax),xtmp(imax));
plot(xtmp(imax),ytmp(imax),'ko')







% % theta1 = 0.7*2*pi;
% theta1 = 0.78*2*pi;
ninterp = 100;
rr = linspace(0,5,ninterp);
xx = rr.*cos(theta1);
yy = rr.*sin(theta1);
if plotit
    plot(xx,yy,'k')
end
% Interpolate along line

dd = interp2(xinterp,yinterp,data,xx,yy);
if plotit
    figure; hold on; box on
    plot(rr,dd)
    xlabel('r [cm]')
    ylabel('\DeltaT')
end

function f = minfun(x)


    dtmp = interp2(xinterp,yinterp,data,xevals+x(1),yevals+x(2));
    f = 1./sum(sum(dtmp));
%     plot(xevals+x(1),yevals+x(2),'g.')    
%     fprintf('x,y,f %e, %e, %e\n',x(1),x(2),f)
end

end



