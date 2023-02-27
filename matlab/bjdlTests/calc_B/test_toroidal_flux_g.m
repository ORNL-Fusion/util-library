clearvars;
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Haskey\g179477.02280';

g = readg_g3d(gfile_name);


psiThis = calc_psi(g,2.1,0)
% psiThis = -0.1;

% nRefines = -2:7;
% for i = 1:length(nRefines)

nRefine = 4;
% nRefine = nRefines(i)
[rPlot,zPlot,psiPlot] = refine_psi(g.r,g.z,2^nRefine,g);
[r2,z2] = ndgrid(rPlot,zPlot);

in = psiPlot < psiThis;
dA = (rPlot(2)-rPlot(1))*(zPlot(2)-zPlot(1));
b = bfield_geq_bicub(g,r2(in),z2(in));
psi_tor = sum(b.bphi)*dA;


if 0
figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14,'fontweight','bold'); grid on; set(gcf,'color','w');
plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),'k-','linewidth',2)
contour(rPlot,zPlot,psiPlot.',psiThis*[1,1],'-','linewidth',1);
contour(rPlot,zPlot,psiPlot.',g.ssibry*[1,1],'k-','linewidth',2);
axis equal; axis tight; ax = axis;
plot(g.rmaxis,g.zmaxis,'bo')
plot(r2,z2,'k.')
plot(r2(in),z2(in),'r.')
end


psiInterp = linspace(g.ssimag,psiThis,100);
psiNInterp = (psiInterp - g.ssimag)/(g.ssibry-g.ssimag);
qInterp = interp1(g.pn,g.qpsi,psiNInterp);
Phi = g.ip_sign*trapz(psiInterp,qInterp)*2*pi


Phi_bry = g.ip_sign*trapz(g.qpsi)*2*pi*(g.ssibry-g.ssimag)/(g.mw-1)

rho = sqrt(Phi/Phi_bry)


Reff = sqrt(abs(2*pi*Phi/(pi*g.bcentr)))
Reff_bry = sqrt(abs(2*pi*Phi_bry/(pi*g.bcentr)))
rhofromReff = Reff/Reff_bry

% q = dPhi/dPsi

%%
function [r2,z2,p2] = refine_psi(r,z,fac,g)
r2 = linspace(r(2),r(end-1),length(r)*fac);
z2 = linspace(z(2),z(end-1),length(z)*fac);
for i = 1:length(z2)
    p2(:,i) = calc_psi(g,r2,z2(i)*ones(size(r2)));
end
end
