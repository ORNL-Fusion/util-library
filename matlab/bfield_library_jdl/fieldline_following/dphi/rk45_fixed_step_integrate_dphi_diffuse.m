function [yout,xout,ierr,i_last_good] = rk45_fixed_step_integrate_dphi_diffuse(y0,x0,dx,nsteps,bfield,nowarn,dmag)
% Same as normal version but applies cross-field kicks
% dmag: Magnetic diffusivity (m^2/m)
if nargin < 6
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
    [dydx,ierr_deriv] = choose_fl_derivs_dphi(x,y,bfield,nowarn);
    if ierr_deriv == 1
        ierr = 1;
        i_last_good = i;
        return;
    end
    [ytmp,ierr_rk4core] = rk4_core_dphi(y,dydx,x,dx,bfield,nowarn);
    if ierr_rk4core == 1
        ierr = 1;
        i_last_good = i;
        return;
    end

    % Step in the current direction
    x = x + dx;
    
    % Diffuse       

    [Bout,ierrB] = bfield_general_rzphi(ytmp(1),ytmp(2),x,bfield,nowarn);
    if ierrB == 1
        ierr = 1;
        i_last_good = i;
        return;
    end        

    dL = sqrt(ytmp(1)^2 + y(1)^2 - 2*ytmp(1)*y(1)*cos(dx) + ytmp(2)^2 + y(2)^2 - 2*ytmp(2)*y(2));

    % B cross z^hat
    perpdir1(1) = Bout.bphi; % r
    perpdir1(2) = -Bout.br; % phi
    perpdir1(3) = 0; % z
    perpdir1 = perpdir1/sqrt(Bout.bphi^2 + Bout.br^2);

    % B cross r^hat
    perpdir2(1) = 0;
    perpdir2(2) = Bout.bz;
    perpdir2(3) = -Bout.bphi;
    perpdir2 = perpdir2/sqrt(Bout.bz^2 + Bout.bphi^2);    

    alpha = 2*pi*rand; % 0 to 2pi
    dca = cos(alpha);
    dsa = sin(alpha);

    delta_x = sqrt(dmag*dL);
    
    ytmp(1)   = ytmp(1)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1));
    x         = x         + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2));
    ytmp(2)   = ytmp(2)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3));

    % Update
    yout(i+1,:) = ytmp;
    xout(i+1) = x;
    y = ytmp;
end
ierr = 0;
i_last_good = nsteps + 1;