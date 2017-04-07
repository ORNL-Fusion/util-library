function [rtarg_map,psi_eval] = map_pt_to_target_mpex(coil,current,ztarg,reval,zeval)
% Maps the point reval,zeval [m] to a radius at the axial position ztarg in [m].
% Interpolate along target radius
% Accepts vector inputs for reval,zeval.
% J.D. Lore
RMAX = min(coil.rwind); % Inner radius of coil
if reval > RMAX
    error('raval > inner radius of coil -- probably does not make sense')
end
if any(reval < 0) || any(zeval < 0)
    error('reval or zeval < 0')
end
if length(reval) ~= length(zeval)
    error('reval and zeval must be same length')
end

ninterp = 500;
zinterp = ztarg*ones(1,ninterp);
rinterp = linspace(0,RMAX,ninterp);
psiinterp = calc_psi_mpex(coil,current,rinterp,zinterp);

% figure; hold on; box on;
% plot(rinterp,psiinterp)

psi_eval = calc_psi_mpex(coil,current,reval,zeval);
rtarg_map = interp1(psiinterp,rinterp,psi_eval);
