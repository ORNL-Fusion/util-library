function [s,ierr,i_last_good]=follow_fieldlines_dphi_gfile_mex(bfield,Rstart,Zstart,phistart,dphi,nsteps,nowarn)
ierr = check_bfield_struct(bfield);
if ierr ~= 0
    error('Cannot follow fieldline, bfield structure not properly set up')
end
if ~strcmp(bfield.type,'gfile')
    error('follow_fieldlines_dphi_gfile_mex only supports bfield.type = ''gfile''')
end
if exist(['calc_psi_mex.',mexext],'file') ~= 3
    error('follow_fieldlines_dphi_gfile_mex requires calc_psi_mex to be compiled and on the MATLAB path')
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

Neq = 2;
Nsys = length(Rstart);
y = zeros(1,Nsys*Neq);
y(1:Neq:Nsys*Neq-1) = Rstart;
y(2:Neq:Nsys*Neq)   = Zstart;
x = phistart;
dx = dphi;

[yout,xout,ierr_rk4,i_last_good] = rk4_fixed_step_integrate_dphi_gfile_mex(y,x,dx,nsteps,bfield,nowarn);

s.r = yout(:,1:Neq:Nsys*Neq-1);
s.z = yout(:,2:Neq:Nsys*Neq);
s.phi = xout;

ierr = ierr_rk4;
