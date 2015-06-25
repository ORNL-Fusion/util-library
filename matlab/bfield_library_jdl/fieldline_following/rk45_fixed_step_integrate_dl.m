function [yout,xout,ierr,i_last_good] = rk45_fixed_step_integrate_dl(y0,x0,dx,nsteps,g,rmp,nowarn,type)
if nargin < 8
    type = 'g3d';
end
if nargin < 7
    nowarn = 0;
end

N = length(y0);

% store initial points
yout = NaN(nsteps+1,N);
xout = NaN(nsteps+1,1);
yout(1,:) = y0;
xout(1) = x0;

y = y0;
x = x0;
for i = 1:nsteps
    if strcmp(type,'g3d')
        [dydx,ierr_deriv] = fl_derivs_dl_g3d(x,y,g,rmp,nowarn);
    elseif strcmp(type,'vmec')
        wout = g;
        [dydx,ierr_deriv] = fl_derivs_dl_vmec(x,y,wout,nowarn);
    elseif strcmp(type,'bspline')
        [dydx,ierr_deriv] = fl_derivs_dl_bspline(x,y,g,rmp,nowarn);
    elseif strcmp(type,'just_coils')
        [dydx,ierr_deriv] = fl_derivs_dl_just_coils(x,y,rmp,nowarn);
    else
        error('bad method')
    end
    if ierr_deriv == 1
        ierr = 1;
        i_last_good = i;
        return;
    end
    [ytmp,ierr_rk4core] = rk4_core_dl(y,dydx,x,dx,g,rmp,nowarn,type);
    if ierr_rk4core == 1
        ierr = 1;
        i_last_good = i;
        return;
    end

    x = x + dx;
    
    yout(i+1,:) = ytmp;
    xout(i+1) = x;
    y = ytmp;
end
ierr = 0;
i_last_good = nsteps + 1;