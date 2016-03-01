function [yout,ierr] = rk4_core(y,dydx,x,dx,bfield,nowarn)
if nargin < 8
    type = 'g3d';
end
if nargin < 7
    nowarn = 0;
end
d1 = dx*dydx;
xtmp = x+dx/2;
ytmp = y+d1/2;
[dydx,ierr_deriv] = choose_fl_derivs(xtmp,ytmp,bfield,nowarn);

% if strcmp(bfield.type,'g3d')
%     [dydx,ierr_deriv] = fl_derivs_dphi_g3d(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'vmec')
%     [dydx,ierr_deriv] = fl_derivs_dphi_vmec(xtmp,ytmp,bfield.g,nowarn);
% elseif strcmp(bfield.type,'bspline')
%     [dydx,ierr_deriv] = fl_derivs_dphi_bspline(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'just_coils')
%     [dydx,ierr_deriv] = fl_derivs_dphi_just_coils(xtmp,ytmp,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'m3dc1')
%     [dydx,ierr_deriv] = fl_derivs_dphi_m3dc1(xtmp,ytmp,bfield.rmp,nowarn);
% else 
%     error('bad method')
% end

if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d2 = dx*dydx;
xtmp = x+dx/2;
ytmp = y+d2/2;
[dydx,ierr_deriv] = choose_fl_derivs(xtmp,ytmp,bfield,nowarn);
% if strcmp(bfield.type,'g3d')
%     [dydx,ierr_deriv] = fl_derivs_dphi_g3d(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'vmec')
%     [dydx,ierr_deriv] = fl_derivs_dphi_vmec(xtmp,ytmp,bfield.g,nowarn);
% elseif strcmp(bfield.type,'bspline')
%     [dydx,ierr_deriv] = fl_derivs_dphi_bspline(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);   
% elseif strcmp(bfield.type,'just_coils')
%     [dydx,ierr_deriv] = fl_derivs_dphi_just_coils(xtmp,ytmp,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'m3dc1')
%     [dydx,ierr_deriv] = fl_derivs_dphi_m3dc1(xtmp,ytmp,bfield.rmp,nowarn);
% end
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d3 = dx*dydx;
xtmp = x+dx;
ytmp = y+d3;
[dydx,ierr_deriv] = choose_fl_derivs(xtmp,ytmp,bfield,nowarn);
% if strcmp(bfield.type,'g3d')
%     [dydx,ierr_deriv] = fl_derivs_dphi_g3d(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'vmec')
%     [dydx,ierr_deriv] = fl_derivs_dphi_vmec(xtmp,ytmp,bfield.g,nowarn);
% elseif strcmp(bfield.type,'bspline')
%     [dydx,ierr_deriv] = fl_derivs_dphi_bspline(xtmp,ytmp,bfield.g,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'just_coils')
%     [dydx,ierr_deriv] = fl_derivs_dphi_just_coils(xtmp,ytmp,bfield.rmp,nowarn);
% elseif strcmp(bfield.type,'m3dc1')
%     [dydx,ierr_deriv] = fl_derivs_dphi_m3dc1(xtmp,ytmp,bfield.rmp,nowarn);
% end
if ierr_deriv == 1
    if ~nowarn
        warning('fl deriv error in rk4_core')
    end
    yout = []; ierr = 1;
    return;
end
d4 = dx*dydx;
yout = y + (d1 + 2*d2 + 2*d3 + d4)./6;
ierr = 0;