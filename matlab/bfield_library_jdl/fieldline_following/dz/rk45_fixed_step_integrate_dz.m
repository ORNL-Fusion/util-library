function [yout,xout,ierr,i_last_good] = rk45_fixed_step_integrate_dz(y0,x0,dx,nsteps,bfield,nowarn)
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
    [dydx,ierr_deriv] = choose_fl_derivs_dz(x,y,bfield,nowarn);
    if ierr_deriv == 1
        ierr = 1;
        i_last_good = i;
        return;
    end
    [ytmp,ierr_rk4core] = rk4_core_dz(y,dydx,x,dx,bfield,nowarn);
    if ierr_rk4core == 1
        ierr = 1;
        i_last_good = i;
        return;
    end

    x = x + dx;
    
    % Here check in vessel?
    if isfield(bfield,'stop_at_vessel')
        if bfield.stop_at_vessel
            if length(bfield.vessel_clip_phi) == 1
                rtest = ytmp(1:2:end);
                ztest = x*ones(size(rtest));
                isin = inpolygon(rtest,ztest,bfield.vessel_clip_r,bfield.vessel_clip_z);
            else
                error('need to write phi vessel interpolation')
            end
            rtest(~isin) = NaN;
            ytmp(1:2:end) = rtest;
            ytmp(2:2:end) = rtest;
        end
    end
    
    yout(i+1,:) = ytmp;
    xout(i+1) = x;
    y = ytmp;
end
ierr = 0;
i_last_good = nsteps + 1;