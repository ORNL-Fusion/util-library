function [rr,dd,radius,angle,x0_final,y0_final] = fit_IR_data(shot,plotit,x0_guess,y0_guess,force)
% [rr,dd,radius,angle,x0_final,y0_final] = fit_IR_data(shot,plotit,x0_guess,y0_guess,force)
% 
% Note: x0,y0_guess in [cm]
%       rr in [cm]
% 
% There are a few ways to run this routine:
%
%  1) Supply x0,y0_guess [cm].
%      - Fitting will be attempted using these values as initial guess.
%  2) Supply x0,y0_guess [cm] and force == 1.
%      - The "guess" is used as the actual center position and only the
%      radius is fit.
%  3) Do not input x0,y0_guess and force.
%      - A window will be opened and you will be prompted to select the
%      approximate center point using the mouse.
%  4) x0,y0_guess == 0 and force is not input.
%      - find_IR_center will be called and used to supply an initial guess.
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
opt_display = 'off';     % Controls optimization output.  'off','final','iter'
ax_size = 7;     % plot axes are ax_size*[-1,1,-1,1]
optimizer = 1;   % 0:fmincon, 1:lsqnonlin w/ LM, 2:lsqnonlin w/ trust-region-reflective
debug = 0;       % display function evaluations
radius_fit_power = 1;  % Mimizing function is weighted by r to this power 
                         % to avoid a fit with a tiny radius in the maximum deltaT region.  
                         % range of 0.7 to 1 seems to work pretty well.

% Read data and create intial cells
fname = find_IR_file(shot);
data = load_IR_file(fname,debug_plots);
cells = create_IR_cells(data,x0_guess,y0_guess);

% Attempt to find center
if manual_select == 0 && x0_guess == 0 && y0_guess == 0
    fprintf('Attempting to automatically find center\n')
    [x0_guess,y0_guess] = find_IR_center(cells,data);
    cells = create_IR_cells(data,x0_guess,y0_guess);
end

% Plot with initial center guess.  If no guess prompt for mouse select.
if plotit
    figure;hold on; box on;
    patch(cells.xcell,cells.ycell,data.IRdata1D,'edgecolor','none')
    colorbar;
    xlabel('X [cm]','fontsize',14)
    ylabel('Y [cm]','fontsize',14)
    title(['\DeltaT, ',num2str(shot)],'fontsize',14)
    axis tight;
    axis(ax_size*[-1,1,-1,1])

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
        title('\DeltaT [C]','fontsize',14)
        axis tight;
        axis(ax_size*[-1,1,-1,1])
    end   
end
fprintf('Using x0_guess = %f, y0_guess = %f [cm]\n',x0_guess,y0_guess)


tevals = linspace(0,2*pi,160); tevals(end) = [];
costt = cos(tevals);
sintt = sin(tevals);

% dx_step = 1;
radius_eval = 1.5;  % initial guess [cm]
lb = [0.1,-10,-10];
ub = [6,10,10];
if force == 1
    x00=radius_eval;
    lb = lb(1);
    ub = ub(1);
else
    x00=[radius_eval,0,0];
end

if optimizer == 0
    fprintf('Using fmincon\n')
    opts=optimoptions(@fmincon,'Display',opt_display); %,'findiffrelstep',dx_step,
    xfinal = fmincon(@maxfun,x00,[],[],[],[],lb,ub,[],opts);
elseif optimizer == 1
    fprintf('Using lsqnonlin Levenberg-Marquardt\n')
    opts=optimoptions(@lsqnonlin,'algorithm','levenberg-marquardt','Display',opt_display);
    xfinal = lsqnonlin(@minfun,x00,[],[],opts);
elseif optimizer == 2
    fprintf('Using lsqnonlin trust-region-reflective\n')
    opts=optimoptions(@lsqnonlin,'Display',opt_display);
    xfinal = lsqnonlin(@minfun,x00,[],[],opts);    
else    
    error('Bad value for optimizer %d',optimizer)
end
radius  = xfinal(1);
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
if optimizer == 0
    fprintf('Found max function value of %f\n',minfun(xfinal))
elseif optimizer > 0
    f = 1/minfun(xfinal);    
    fprintf('Found max function value of %f\n',f)
end

% Use fit to evaluate
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
    axis(ax_size*[-1,1,-1,1])
end
theta = linspace(0,2*pi,100);
if plotit
    plot(radius*cos(theta),radius*sin(theta),'m')
end

% 
xmean = mean(cells.xcell);
ymean = mean(cells.ycell);
xeval = 1.2*radius*cos(theta);
yeval = 1.2*radius*sin(theta);
isin = inpolygon(xmean,ymean,xeval,yeval);
[dmax,imax] = max(data.IRdata1D(isin));
xtmp =  xmean(isin);
ytmp =  ymean(isin);
angle = atan2(ytmp(imax),xtmp(imax));
if plotit
    plot(xtmp(imax),ytmp(imax),'ko')
end
ninterp = 100;
rr = linspace(0,5,ninterp);
xx = rr.*cos(angle);
yy = rr.*sin(angle);
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

    function f = maxfun(x)
        f = -minfun(x);
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
        
        if optimizer == 0
            f  = sum(dtmp)*x(1)^radius_fit_power;
        elseif optimizer > 0
            f = 1/(x(1)^radius_fit_power*sum(dtmp));      
        end

            if debug && length(x) == 3 && length(f) == 1
                fprintf('x: %f %f %f, y: %f\n',x,f)
            elseif debug && length(x) == 1 && length(f) == 1
                fprintf('x: %f, y: %f\n',x,f)
            end
        %     f = [1,1./x(1)]./sum(sum(dtmp));
        
        %     f = 1./dtmp;

        if debug_plots >= 2
            plot(xevals-xshift,yevals-yshift,'g.')
        end
    end




end




