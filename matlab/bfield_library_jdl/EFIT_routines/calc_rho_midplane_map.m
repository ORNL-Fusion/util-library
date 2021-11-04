function [rho,R,Rsep] = calc_rho_midplane_map(g,R1,Z1)
% Give R,Z or just R = psiN
% rho = R - Rsep

if nargin == 2
    psiN1 = R1;
end

Ninterp = 100;
Rinterp = linspace(g.rmaxis,g.r(end-1),Ninterp);
Zinterp = g.zmaxis*ones(1,Ninterp);
[psiN_mid,psi_mid] = calc_psiN(g,Rinterp,Zinterp);

% plot_gfile(g)
% plot(Rinterp,Zinterp,'cx')
% plot(R1,Z1,'cx','markersize',12,'linewidth',3)


if nargin == 3
    [psiN1,psi1] = calc_psiN(g,R1,Z1);
end
R_midplane_map = interp1(psiN_mid,Rinterp,psiN1);   

Rsep = interp1(psiN_mid,Rinterp,1);   
rho_midplane_map = R_midplane_map - Rsep;

rho = rho_midplane_map;
R = R_midplane_map;

% figure; hold on;
% plot(Rinterp,psiN_mid)
% plot(R_midplane_map,psiN1,'ro')
