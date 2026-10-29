function [yout,xout,ierr,i_last_good] = rk4_fixed_step_integrate_dphi_gfile_mex(y0,x0,dx,nsteps,bfield,nowarn)
if nargin < 6
    nowarn = 0;
end

N = length(y0);

yout = NaN(nsteps+1,N);
xout = NaN(nsteps+1,1);
yout(1,:) = y0;
xout(1) = x0;

y = y0;
x = x0;
for i = 1:nsteps
    [dydx,ierr_deriv] = fl_derivs_dphi_gfile_mex(y,bfield,nowarn);
    if ierr_deriv == 1
        ierr = 1;
        i_last_good = i;
        return;
    end

    d1 = dx*dydx;

    ytmp = y + d1/2;
    [dydx,ierr_deriv] = fl_derivs_dphi_gfile_mex(ytmp,bfield,nowarn);
    if ierr_deriv == 1
        if ~nowarn
            warning('fl deriv error in rk4_fixed_step_integrate_dphi_gfile_mex')
        end
        ierr = 1;
        i_last_good = i;
        return;
    end
    d2 = dx*dydx;

    ytmp = y + d2/2;
    [dydx,ierr_deriv] = fl_derivs_dphi_gfile_mex(ytmp,bfield,nowarn);
    if ierr_deriv == 1
        if ~nowarn
            warning('fl deriv error in rk4_fixed_step_integrate_dphi_gfile_mex')
        end
        ierr = 1;
        i_last_good = i;
        return;
    end
    d3 = dx*dydx;

    ytmp = y + d3;
    [dydx,ierr_deriv] = fl_derivs_dphi_gfile_mex(ytmp,bfield,nowarn);
    if ierr_deriv == 1
        if ~nowarn
            warning('fl deriv error in rk4_fixed_step_integrate_dphi_gfile_mex')
        end
        ierr = 1;
        i_last_good = i;
        return;
    end
    d4 = dx*dydx;

    ytmp = y + (d1 + 2*d2 + 2*d3 + d4)./6;
    x = x + dx;

    yout(i+1,:) = ytmp;
    xout(i+1) = x;
    y = ytmp;
end

ierr = 0;
i_last_good = nsteps + 1;
