function [s,ierr,i_last_good]=follow_fieldlines_rzphi_dl(g,rmp,Rstart,Zstart,phistart,dl,nsteps,nowarn,type)
% if type == vmec, then g = wout, rmp is not used
% Phi in radians, RZ in meters
if nargin < 9
    type = 'g3d';
end
if nargin < 8
    nowarn = 0;
end

Neq = 3;                % Two equations for each ode system 
Nsys = length(Rstart);  % Number of simultaneous systems to be solved
y(1:Neq:Nsys*Neq-2) = Rstart;
y(2:Neq:Nsys*Neq-1) = phistart;
y(3:Neq:Nsys*Neq)   = Zstart;
x = 0;
dx = dl;

[yout,xout,ierr_rk45,i_last_good] = rk45_fixed_step_integrate_dl(y,x,dx,nsteps,g,rmp,nowarn,type);

s.r   = yout(:,1:Neq:Nsys*Neq-2);
s.phi = yout(:,2:Neq:Nsys*Neq-1);
s.z   = yout(:,3:Neq:Nsys*Neq);
s.l  = xout;

ierr = ierr_rk45;







