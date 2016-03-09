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
debug_plots = 0; % Turn on for plots. Can be level [0,1,2,3]
verbose = 0;     % Controls lsqnonlin output.  Can be level [0,1,2]


fname = find_IR_file(shot);
data = load_IR_file(fname,debug_plots);



cells = create_IR_cells(data,x0_guess,y0_guess);
if manual_select == 0 && x0_guess == 0 && y0_guess == 0
    fprintf('Attempting to automatically find center\n')
    [x0_guess,y0_guess] = find_IR_center(cells,data);
    cells = create_IR_cells(data,x0_guess,y0_guess);
end

tolx = 1e-8;
tolfun = 1e-8;

% Plot with initial center guess.  If no guess ask for click
if plotit
    figure;hold on; box on;
    patch(cells.xcell,cells.ycell,data.IRdata1D,'edgecolor','none')
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

        cells = create_IR_cells(data,x0_guess,y0_guess);
        clf;hold on; box on;
        patch(cells.xcell,cells.ycell,data.IRdata1D,'edgecolor','none')
        colorbar;
        xlabel('X [cm]','fontsize',14)
        ylabel('Y [cm]','fontsize',14)
        title('\DeltaT','fontsize',14)
        axis tight;
        axis([-15,15,-15,15])
    end   
end

% Do the fit

if verbose == 0
    qval = 'off';
elseif verbose == 1
    qval = 'final';
elseif verbose == 2
    qval = 'iter';
else
    qval = 'final';
end
opts=optimoptions(@lsqnonlin,'TolFun',tolfun,'TolX',tolx,'Display',qval);
min_method = 0;  % 0 = LM, 1 = trust
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
    A = 1;
    b = 10;
else
    x00=[radius_eval,0,0];
    A = [1,0,0];
    b = 10;
end

xfinal = fmincon(@minfun,x00,A,b);

% xfinal=lsqnonlin(@minfun,x00,lb,ub,opts);
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

cells = create_IR_cells(data,x0_final,y0_final);

if debug_plots >= 2 && plotit
    figure;hold on; box on;
elseif plotit 
    clf; hold on; box on;    
end
if plotit    
    patch(cells.xcell,cells.ycell,data.IRdata1D,'edgecolor','none')
    plot([0,0],[min(min(cells.ycell)),max(max(cells.ycell))],'k-')
    plot([min(min(cells.xcell)),max(max(cells.xcell))],[0,0],'k-')
    colorbar;
    xlabel('X [cm]','fontsize',14)
    title(['\DeltaT, ',num2str(shot)],'fontsize',14)
    ylabel('Y [cm]','fontsize',14)
    axis tight;
%     axis(15*[-1,1,-1,1])
    axis(5*[-1,1,-1,1])
end

theta = linspace(0,2*pi,100);
if plotit
    plot(radius*cos(theta),radius*sin(theta),'m')
end

xmean = mean(cells.xcell);
ymean = mean(cells.ycell);
xeval = 1.2*radius*cos(theta);
yeval = 1.2*radius*sin(theta);
isin = inpolygon(xmean,ymean,xeval,yeval);
[dmax,imax] = max(data.IRdata1D(isin));
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
dd = interp2(cells.xinterp1D,cells.yinterp1D,data.IRdata2D,xx,yy);
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
    dtmp = interp2(cells.xinterp1D,cells.yinterp1D,data.IRdata2D,xevals+xshift,yevals+yshift);
%     f = -sum(dtmp);
%     if length(x) == 3 && length(f) == 1
%         fprintf('x: %f %f %f, y: %f\n',x,f)
%     else
%         fprintf('x: %f, y: %f\n',x,f)
%     end
%     f = [1,1./x(1)]./sum(sum(dtmp));
        f  = -sum(dtmp);
%     f = 1./dtmp;
%     f = -sum(abs(dtmp));
    if debug_plots >= 2
        plot(xevals-xshift,yevals-yshift,'g.')
    end
end




end




