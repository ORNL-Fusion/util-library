function phi_tor = calculate_toroidal_flux(bfield,R_surf,Z_surf,phi_surf)
%-------------------------------------------------------------------

tic;
% R_surf = poinc.Rpoinc;
% Z_surf = poinc.Zpoinc;
% phi = poinc.phi;
phi = phi_surf;

Rax = mean(R_surf);
Zax = mean(Z_surf);
disp(['Guess of axis from geometric average of lcfs coords: Rax=' num2str(Rax) ' Zax=' num2str(Zax)])


% Set integration limits
dx = 1e-3;
Rmin = min(R_surf) - dx;
Rmax = max(R_surf) + dx;
Zmin = min(Z_surf) - dx;
Zmax = max(Z_surf) + dx;

theta_surf = atan2(Z_surf-Zax,R_surf-Rax);
r_surf = sqrt((R_surf-Rax).^2 + (Z_surf-Zax).^2);
[theta_surf,Isort]=sort(theta_surf);                       %sort by increasing theta
r_surf = r_surf(Isort);
Z_surf = Z_surf(Isort);
R_surf = R_surf(Isort);


nx = 100;
Rg = linspace(Rmin,Rmax,nx);
Zg = linspace(Zmin,Zmax,nx);
dR = Rg(2)-Rg(1);
dZ = Zg(2)-Zg(1);
flux2 = 0;
nowarn = 0;

figure; hold on; box on;
plot(R_surf,Z_surf,'k.-')
for rind = 1:length(Rg)
    for zind = 1:length(Zg)
        ikeep = inpolygon(Rg(rind),Zg(zind),R_surf,Z_surf);
        if ikeep
            [Bout,ierr] = bfield_general_rzphi(Rg(rind),Zg(zind),phi,bfield,nowarn);
            flux2 = flux2 + Bout.bphi*dR*dZ;
            plot(Rg(rind),Zg(zind),'r.','MarkerSize',0.5)            
        else
            plot(Rg(rind),Zg(zind),'b.','MarkerSize',0.5)
        end
    end
end
disp(['Flux estimate: ' num2str(flux2)])
fprintf('Flux calc took %f seconds.\n',toc); tic;
phi_tor = flux2;