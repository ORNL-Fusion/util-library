function [s,ierr,i_last_good]=follow_fieldlines_rzphi(g,rmp,Rstart,Zstart,phistart,dphi,nsteps,nowarn,type)
% if type == vmec, then g = wout, rmp is not used
% Phi in radians, RZ in meters
if nargin < 9
    type = 'g3d';
end
if nargin < 8
    nowarn = 0;
end
if length(phistart) > 1
    error('wtf')
end

Neq = 2;  % Two equations for each ode system 
Nsys = length(Rstart);  % Number of simultaneous systems to be solved
y(1:Neq:Nsys*Neq-1) = Rstart;
y(2:Neq:Nsys*Neq)   = Zstart;
x = phistart;
dx = dphi;

[yout,xout,ierr_rk45,i_last_good] = rk45_fixed_step_integrate(y,x,dx,nsteps,g,rmp,nowarn,type);

s.r = yout(:,1:Neq:Nsys*Neq-1);
s.z = yout(:,2:Neq:Nsys*Neq);
s.phi = xout;

ierr = ierr_rk45;







