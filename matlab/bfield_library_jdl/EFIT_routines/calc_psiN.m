function [psiN,psi,ierr] = calc_psiN(g,R,Z,quiet)
%[psiN,psi] = calc_psiN(g,R,Z,quiet)
% R, Z can be vectors
% quiet optional

if nargin < 4
    quiet = 0;
end
[psi,ierr]= calc_psi(g,R,Z,quiet);         % this applies g.ip_sign
psiN = (g.ip_sign*psi-g.ssimag)/(g.ssibry-g.ssimag); % so have to apply it again because ssi quantities will be flipped (if ip_sign = -1)