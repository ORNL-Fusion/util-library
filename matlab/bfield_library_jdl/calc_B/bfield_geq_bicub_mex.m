function [Bout,ierr] = bfield_geq_bicub_mex(g,R,Z,nowarn)
% Sidecar version of bfield_geq_bicub using calc_psi_mex.

if nargin < 4
    nowarn = 0;
end

if iscolumn(R)
    R = R';
    Z = Z';
end

toroidal_off_grid = 0;
if isfield(g,'toroidal_off_grid')
    toroidal_off_grid = g.toroidal_off_grid;
end

[psi,ierr,dpsidr,dpsidz] = calc_psi_mex(g,R,Z,nowarn);

R1_inv = 1./R;
br1 = -dpsidz.*R1_inv;
bz1 =  dpsidr.*R1_inv;

psiN = (g.ip_sign*psi-g.ssimag)/(g.ssibry-g.ssimag);

% Use the tabulated current function on the core branch and vacuum field outside.
bt1 = g.bcentr*g.rzero.*R1_inv;
ic = find(psiN <= 1.0);
if ~isempty(ic)
    fpol = ppval(g.fpol_pp,psiN(ic));
    bt1(ic) = fpol.*R1_inv(ic);
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
