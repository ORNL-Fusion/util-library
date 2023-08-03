function [Bout,ierr] = bfield_equ_bicub(Equ,R,Z,nowarn)
% [Bout,ierr] = bfield_equ_bicub(g,R1,Z1,nowarn)
% Bout is nan when R or Z is off grid
% Uses format from read_equ

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
if isfield(Equ,'toroidal_off_grid')
    toroidal_off_grid = Equ.toroidal_off_grid;
end        

if any(isnan(R)) || any(isnan(Z))
    error('R or Z is NaN')
end

[psi,ierr,dpsidr,dpsidz] = calc_psi(Equ,R,Z,nowarn);
    
R1_inv = 1./R;
br1 = -dpsidz.*R1_inv;
bz1 =  dpsidr.*R1_inv;

% toroidal field
bt1 = Equ.btf*Equ.rtf.*R1_inv;

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


