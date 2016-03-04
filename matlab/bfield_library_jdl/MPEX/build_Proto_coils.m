function [coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose)
if nargin < 5
    verbose = 1;
end
if nargin < 4
    error('must specify inputs')
end
debug_plots = 0;

[nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils(helicon_current,current_A,current_B,config,verbose);

ntheta_per_wind = 100;

ncoils = length(nturns);
for i=1:ncoils
    [coil0,current0] = build_circular_coil(rr1(i),rr2(i),z0(i),cl(i),nturns(i),nlayers(i),cur(i),ntheta_per_wind);
    if i > 1
        coil = [coil;coil0];
        current = [current;current0];
    else
        coil = coil0;
        current = current0;
    end  
end

if debug_plots
    figure; hold on; box on;
    % plot3(coil(:,1),coil(:,2),coil(:,3))
    for i = 1:length(current)-1
        if current(i) < 1e-8
            plot3(coil(i:i+1,1),coil(i:i+1,2),coil(i:i+1,3),'r')
        else
            plot3(coil(i:i+1,1),coil(i:i+1,2),coil(i:i+1,3),'b')
        end
    end    
end