function [psiN,psi] = calc_psiN(g,R1,Z1,quiet)
%[psiN,psi] = calc_psiN(g,R1,Z1,quiet)
% R1, Z1 can be vectors
% quiet optional
if nargin < 4
    quiet = 0;
end
psi = g.ip_sign*get_psi_bicub(g,R1,Z1,quiet);
psiN = (psi-g.ssimag)/(g.ssibry-g.ssimag);