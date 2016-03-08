function [rr,dd,radius,angle,x0_final,y0_final] = plot_IR_data_raw(shot,plotit,x0_guess,y0_guess,force)
% [rr,dd,radius,angle,x0_final,y0_final] = plot_IR_data_raw(shot,plotit,x0_guess,y0_guess,force)
% x0_guess, y0_guess [cm]. Guess of center of plasma from the middle of the
% frame.  For the 74xx shots I've looked at [0,-2.5] works well.
%   
% if force == 1 the guess is not changed
%  
if nargin < 2
    plotit = 1;
end
if nargin < 3
    plotit = 1;
    manual_select = 1;
    x0_guess = 0;
    y0_guess = 0;
else
    manual_select = 0;
end
if nargin < 5
    force = 0;
end

% Settings
px_per_cm = 12.146;
debug_plots = 0; % Turn on for plots. Can be level [0,1,2,3]
verbose = 0;     % Controls lsqnonlin output.  Can be level [0,1,2]
tolx = 1e-8;
tolfun = 1e-8;

% Search data path for file names including string of shot number
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

% plot raw data
if debug_plots >= 3
    figure; hold on; box on;
    h = pcolor(d.Frame);
    colorbar;
    set(h,'edgecolor','none')
    title('Raw data','fontsize',14)
    xlabel('X [px]','fontsize',14)
    ylabel('Y [px]','fontsize',14)
    set(gca,'fontsize',14)
end

% Convert to cm and flip
data = flipud(d.Frame);
nw = size(d.Frame,2);
nh = size(d.Frame,1);
dx = nw/px_per_cm;
dy = nh/px_per_cm;
icount = 1;
dcell = zeros(1,nh*nw);
for i = 1:nw
    for j = 1:nh
        dcell(icount) = data(j,i);
        icount = icount + 1;
    end
end
[x,y,xinterp,yinterp,xcell,ycell,xmesh,ymesh] = create_cells(dx,dy,nh,nw,x0_guess,y0_guess);

% Plot with initial center guess.  If no guess ask for click
if plotit
    figure;hold on; box on;
    patch(xcell,ycell,dcell,'edgecolor','none')
    colorbar;
    xlabel('X [cm]','fontsize',14)
    ylabel('Y [cm]','fontsize',14)
    title(['\DeltaT, ',num2str(shot)],'fontsize',14)
    axis tight;
    axis([-15,15,-15,15])

    if manual_select
        title('Click on approximate center position')
        fprintf('Waiting for mouse input\n')
        [x0_guess,y0_guess] = ginput(1);
        fprintf('Click gave x0_guess = %f, y0_guess = %f\n',x0_guess,y0_guess)

        [x,y,xinterp,yinterp,xcell,ycell,xmesh,ymesh] = create_cells(dx,dy,nh,nw,x0_guess,y0_guess);
        clf;hold on; box on;
        patch(xcell,ycell,dcell,'edgecolor','none')
        colorbar;
        xlabel('X [cm]','fontsize',14)
        ylabel('Y [cm]','fontsize',14)
        title('\DeltaT','fontsize',14)
        axis tight;
        axis([-15,15,-15,15])
    end
    % x = linspace(0,nw/px_per_cm);
    % y = linspace(0,nh)/px_per_cm;
    % h= pcolor(x,y,flipud(d.Frame));
    % set(h,'edgecolor','none')    
end

% Perform fit for center and radius through maximum
if verbose == 0
    qval = 'off';
elseif verbose == 1
    qval = 'final';
elseif verbose == 2
    qval = 'iter';
else
    qval = 'final';
end
opts=optimoptions(@lsqnonlin,'TolFun',tolfun,'TolX',tolx,'Display',qval,'FinDiffType','central','diffminchange',.1); %,'TypicalX',20000000*[1,1]);

% opts=optimoptions(@lsqnonlin);
min_method = 0;  % 0 = LM, 1 = trust


%  try to guess center


if min_method == 0
    opts.Algorithm = 'levenberg-marquardt';
    lb = [];
    ub = [];
elseif min_method == 1
    opts.Algorithm = 'trust-region-reflective';
    lb = [-10,-10];
    ub = [10,10];
end

radius_eval = 9;  %cm
x00=[0,0];
xfinal=lsqnonlin(@minfun2,x00,lb,ub,opts);
x0_final = xfinal(1) + x0_guess;
y0_final = xfinal(2) + y0_guess;


if plotit 
    plot(x0_final,y0_final,'mx')
    plot(x00(1),x00(2),'mo')
end
theta2 = linspace(0,2*pi,100);
plot(radius_eval*cos(theta2)-x0_final,radius_eval*sin(theta2)-y0_final,'m-')
fprintf('Found x0,y0 of %f,%f\n',x0_final,y0_final)

asfadsfadsfasdf




% Do the fit

min_method = 1;  % 0 = LM, 1 = trust
if min_method == 0
    opts.Algorithm = 'levenberg-marquardt';
    lb = [];
    ub = [];
elseif min_method == 1
    opts.Algorithm = 'trust-region-reflective';
    lb = [0.1,-10,-10];
    ub = [5,10,10];
end

radius_eval = 1.5;  %cm
tevals = linspace(0,2*pi,80); tevals(end) = [];
costt = cos(tevals);
sintt = sin(tevals);

if force == 1
    x00=[radius_eval];
    if ~isempty(lb)
        lb = lb(1);
        ub = ub(1);
    end
else
    x00=[radius_eval,0,0];
end
xfinal=lsqnonlin(@minfun,x00,lb,ub,opts);
radius   = xfinal(1);
if force== 1
    x0_final = x0_guess;
    y0_final = y0_guess;
else
    x0_final = xfinal(2) + x0_guess;
    y0_final = xfinal(3) + y0_guess;
end

if plotit 
    plot(x0_final,y0_final,'mx')
    if force ~= 1
        plot(x00(2),x00(3),'mo')
    end
end

fprintf('Found x0,y0,radius of %f,%f,%f\n',x0_final,y0_final,radius)

[x,y,xinterp,yinterp,xcell,ycell] = create_cells(dx,dy,nh,nw,x0_final,y0_final);

if debug_plots >= 2 && plotit
    figure;hold on; box on;
elseif plotit 
    clf; hold on; box on;    
end
if plotit    
    patch(xcell,ycell,dcell,'edgecolor','none')
    plot([0,0],[min(min(ycell)),max(max(ycell))],'k-')
    plot([min(min(xcell)),max(max(xcell))],[0,0],'k-')
    colorbar;
    xlabel('X [cm]','fontsize',14)
    title(['\DeltaT, ',num2str(shot)],'fontsize',14)
    ylabel('Y [cm]','fontsize',14)
    axis tight;
    axis([-15,15,-15,15])
end

theta = linspace(0,2*pi,100);
if plotit
    plot(radius*cos(theta),radius*sin(theta),'m')
end

xmean = mean(xcell);
ymean = mean(ycell);
xeval = 1.2*radius*cos(theta);
yeval = 1.2*radius*sin(theta);
isin = inpolygon(xmean,ymean,xeval,yeval);
[dmax,imax] = max(dcell(isin));
xtmp =  xmean(isin);
ytmp =  ymean(isin);
theta1 = atan2(ytmp(imax),xtmp(imax));
if plotit
    plot(xtmp(imax),ytmp(imax),'ko')
end
angle = theta1;
ninterp = 100;
rr = linspace(0,5,ninterp);
xx = rr.*cos(theta1);
yy = rr.*sin(theta1);
if plotit
    plot(xx,yy,'k')
end

% Interpolate along line
dd = interp2(xinterp,yinterp,data,xx,yy);
if plotit >= 2
    figure; hold on; box on
    plot(rr,dd)
    ylim = get(gca,'ylim');
    plot(radius*[1,1],ylim,'k')
    title(num2str(shot),'fontsize',14)
    xlabel('r [cm]','fontsize',14)
    ylabel('\DeltaT','fontsize',14)
    set(gca,'fontsize',14)    
end

function f = minfun(x)
    xevals = x(1).*costt;
    yevals = x(1).*sintt;
    if length(x) == 3
        xshift = x(2);
        yshift = x(3);
    else
        xshift = 0;
        yshift = 0;
    end
    dtmp = interp2(xinterp,yinterp,data,xevals+xshift,yevals+yshift);
%     f = [1,1./x(1)]./sum(sum(dtmp));
    f = 1./dtmp;
    if debug_plots >= 2
        plot(xevals-xshift,yevals-yshift,'g.')
    end
end

function f = maxfun(xshift,yshift)
    rmesh = sqrt((xmesh-xshift).^2 + (ymesh-yshift).^2);
    f = sum(sum(data(rmesh<radius_eval)));
    fprintf('x: %f %f, y: %f\n',x,f)
    if debug_plots >= 2
        plot(xevals-xshift,yevals-yshift,'g.')
    end
end


end

function [x,y,xinterp,yinterp,xcell,ycell,xmesh,ymesh] = create_cells(dx,dy,nh,nw,xshift,yshift)
x = linspace(-dx/2,dx/2,nw+1);
y = linspace(-dy/2,dy/2,nh+1);
x = x - xshift;
y = y - yshift;

icount = 1;
xcell = zeros(4,nh*nw);
ycell = zeros(4,nh*nw);
for i = 1:nw
    for j = 1:nh
        xcell(:,icount) = [x(i),x(i+1),x(i+1),x(i)];
        ycell(:,icount) = [y(j),y(j),y(j+1),y(j+1)];
        icount = icount + 1;
    end
end

xinterp = x(1:end-1) + diff(x)/2;
yinterp = y(1:end-1) + diff(y)/2;
[xmesh,ymesh] = meshgrid(xinterp,yinterp);
end


