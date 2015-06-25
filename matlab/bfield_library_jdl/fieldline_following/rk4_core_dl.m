function [yout,ierr] = rk4_core_dl(y,dydx,x,dx,g,rmp,nowarn,type)
if nargin < 8
    type = 'g3d';
end
if nargin < 7
    nowarn = 0;
end
d1 = dx*dydx;
if strcmp(type,'g3d')
    [dydx,ierr_deriv] = fl_derivs_dl_g3d(x+dx/2,y+d1/2,g,rmp,nowarn);
elseif strcmp(type,'vmec')
    wout = g;
    [dydx,ierr_deriv] = fl_derivs_dl_vmec(x+dx/2,y+d1/2,wout,nowarn);
elseif strcmp(type,'bspline')
    [dydx,ierr_deriv] = fl_derivs_dl_bspline(x+dx/2,y+d1/2,g,rmp,nowarn);
elseif strcmp(type,'just_coils')
    [dydx,ierr_deriv] = fl_derivs_dl_just_coils(x+dx/2,y+d1/2,rmp,nowarn);
else 
    error('bad method')
end
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d2 = dx*dydx;
if strcmp(type,'g3d')
    [dydx,ierr_deriv] = fl_derivs_dl_g3d(x+dx/2,y+d2/2,g,rmp,nowarn);
elseif strcmp(type,'vmec')
    wout = g;
    [dydx,ierr_deriv] = fl_derivs_dl_vmec(x+dx/2,y+d2/2,wout,nowarn);
elseif strcmp(type,'bspline')
    [dydx,ierr_deriv] = fl_derivs_dl_bspline(x+dx/2,y+d2/2,g,rmp,nowarn);   
elseif strcmp(type,'just_coils')
    [dydx,ierr_deriv] = fl_derivs_dl_just_coils(x+dx/2,y+d2/2,rmp,nowarn);    
end
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d3 = dx*dydx;
if strcmp(type,'g3d')
    [dydx,ierr_deriv] = fl_derivs_dl_g3d(x+dx,y+d3,g,rmp,nowarn);
elseif strcmp(type,'vmec')
    wout = g;
    [dydx,ierr_deriv] = fl_derivs_dl_vmec(x+dx,y+d3,wout,nowarn);
elseif strcmp(type,'bspline')
    [dydx,ierr_deriv] = fl_derivs_dl_bspline(x+dx,y+d3,g,rmp,nowarn);    
elseif strcmp(type,'just_coils')
    [dydx,ierr_deriv] = fl_derivs_dl_just_coils(x+dx,y+d3,rmp,nowarn);    
end
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d4 = dx*dydx;
yout = y + (d1 + 2*d2 + 2*d3 +d4)./6;
ierr = 0;