function [s,ierr,i_last_good]=follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps,nowarn)
% [s,ierr,i_last_good]=follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps,nowarn)
% Rstart, Zstart in meters (can be arrays)
% phistart in radians, must be scalar
% dphi in radians
% nowarn == 1 supresses warnings from fl_derivs, and rk45
% bfield should be a struct with fields 'type' and others corresponding to the type
%   type = 'gfile'
%       g = gfile info
ierr = check_bfield_struct(bfield);
if ierr ~= 0
    error('Cannot follow fieldline, bfield structure not properly set up')
end
if nargin < 7
    nowarn = 0;
end
if length(phistart) > 1
    error('phistart must be a scalar')
end
if nsteps <= 0
    error('Nsteps must be positive!')
end

Neq = 2;                % Number of equations for each ode system 
Nsys = length(Rstart);  % Number of simultaneous systems to be solved
y(1:Neq:Nsys*Neq-1) = Rstart;
y(2:Neq:Nsys*Neq)   = Zstart;
x = phistart;
dx = dphi;

[yout,xout,ierr_rk45,i_last_good] = rk45_fixed_step_integrate_dphi(y,x,dx,nsteps,bfield,nowarn);

s.r = yout(:,1:Neq:Nsys*Neq-1);
s.z = yout(:,2:Neq:Nsys*Neq);
s.phi = xout;

ierr = ierr_rk45;







