function [Bout,ierr] = bfield_geq_bicub(g,R,Z,nowarn)
% [Bout,ierr] = bfield_geq_bicub(g,R1,Z1,nowarn)
% Bout is nan when R or Z is off grid
% g.ip_sign = -sign(Ip) is applied here
% JDL
%

if nargin < 4
    nowarn = 0;
end

if iscolumn(R)
    R = R';
    Z = Z';
end

% This is a switch to return B = Bt = 1 when the evaluation point(s)
% are off the grid. 
toroidal_off_grid = 0;
if isfield(g,'toroidal_off_grid')
    toroidal_off_grid = g.toroidal_off_grid;
end        

% if any(isnan(R)) || any(isnan(Z))
%     % error('R or Z is NaN')
%     if ~nowarn
%         fprintf('Warning: R or Z is NaN in bfield_geq_bicub\n');        
%     end
% 
% end

[psi,ierr,dpsidr,dpsidz] = calc_psi(g,R,Z,nowarn);
    
R1_inv = 1./R;
br1 = -dpsidz.*R1_inv;
bz1 =  dpsidr.*R1_inv;

% psiN = (psi - g.ssimag)./(g.ssibry-g.ssimag);
psiN = (g.ip_sign*psi-g.ssimag)/(g.ssibry-g.ssimag); 

if any(psiN < 0) 
    error('PsiN out of bounds for fpol interpolation')
end

fpol = polyval(g.fpol_coeffs,psiN);

% toroidal field
bt1 = g.bcentr*g.rzero.*R1_inv;
ic = find(psiN <= 1.0);
if ~isempty(ic)
    bt1(ic) = fpol(ic).*R1_inv(ic);
end

Bout.br = br1;
Bout.bz = bz1;
Bout.bphi = bt1;

Bout.br(ierr) = nan;
Bout.bz(ierr) = nan;
Bout.bphi(ierr) = nan;

if any(ierr)
    if toroidal_off_grid
        if ~nowarn
            warning('Point(s) off grid --> returning toroidal field = 1 there')
        end
        Bout.br(ierr) = 0;
        Bout.bz(ierr) = 0;
        Bout.bphi(ierr) = 1;
        ierr(ierr) = false; 
    end    
end


