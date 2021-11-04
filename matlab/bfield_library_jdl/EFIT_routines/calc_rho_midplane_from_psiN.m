function [rho,R,Rsep] = calc_rho_midplane_from_psiN(g,psiN)

% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\LDRD 2021\d3d data\g175060.02500';
% 
% g = readg_g3d(gfile_name);

% psiN = [0.99,1.01];

nInterp = 100;
rInterp = linspace(g.rmaxis,g.r(end-1),nInterp);
zInterp = g.zmaxis*ones(size(rInterp));

psiNInterp = calc_psiN(g,rInterp,zInterp);

if psiN > max(psiNInterp) | psiN < min(psiNInterp)
    error('PsiN value out of range %f',psiN)
end
R = interp1(psiNInterp,rInterp,psiN);
Rsep = interp1(psiNInterp,rInterp,1);

rho = R - Rsep;

