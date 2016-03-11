function [yout,ierr] = rk4_core_dphi(y,dydx,x,dx,bfield,nowarn)
if nargin < 6
    nowarn = 0;
end
d1 = dx*dydx;
xtmp = x+dx/2;
ytmp = y+d1/2;
[dydx,ierr_deriv] = choose_fl_derivs_dphi(xtmp,ytmp,bfield,nowarn);
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core_dphi')
    end
    yout = []; ierr = 1;
    return;
end
d2 = dx*dydx;
xtmp = x+dx/2;
ytmp = y+d2/2;
[dydx,ierr_deriv] = choose_fl_derivs_dphi(xtmp,ytmp,bfield,nowarn);
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core_dphi')
    end
    yout = []; ierr = 1;
    return;
end
d3 = dx*dydx;
xtmp = x+dx;
ytmp = y+d3;
[dydx,ierr_deriv] = choose_fl_derivs_dphi(xtmp,ytmp,bfield,nowarn);

if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core_dphi')
    end
    yout = []; ierr = 1;
    return;
end
d4 = dx*dydx;
yout = y + (d1 + 2*d2 + 2*d3 + d4)./6;
ierr = 0;