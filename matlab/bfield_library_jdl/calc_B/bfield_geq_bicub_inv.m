function [Bout,ierr] = bfield_geq_bicub_inv(g,R1,Z1,nowarn)
%[Bout,ierr] = bfield_geq_bicub_inv(g,R1,Z1,nowarn)
if nargin < 4
    nowarn = 0;
end

if iscolumn(R1)
    R1 = R1.';
    Z1 = Z1.';
end

% This is a switch to return B = Bt = 1 when the evaluation point(s)
% are off the grid. 
toroidal_off_grid = 0;
if isfield(g,'toroidal_off_grid')
    toroidal_off_grid = g.toroidal_off_grid;
end        

ir = floor( (R1 - g.r(1))/g.dR ) + 1;
iz = floor( (Z1 - g.z(1))/g.dZ ) + 1;

if any(isnan(R1)) | any(isnan(Z1))
    error('R1 or Z1 is NaN')
end

% check for points off grid  (account for derivatives at edge)
iroff = find(ir <= 2 | ir >= g.mw - 1);
izoff = find(iz <= 2 | iz >= g.mh - 1);

mask = zeros(size(R1));
mask(iroff) = 1;
mask(izoff) = 1;
ir(mask==1) = 3;  % just a value that will not error
iz(mask==1) = 3;

dr = (R1 - g.r(ir))/g.dR;
dz = (Z1 - g.z(iz))/g.dZ;

index = ir + (g.mw-1)*(iz-1);
c = g.bicub_coeffs_inv;

c00 = c.c00(index);
c10 = c.c10(index);
c20 = c.c20(index);
c30 = c.c30(index);

c01 = c.c01(index);
c11 = c.c11(index);
c21 = c.c21(index);
c31 = c.c31(index);

c02 = c.c02(index);
c12 = c.c12(index);
c22 = c.c22(index);
c32 = c.c32(index);

c03 = c.c03(index);
c13 = c.c13(index);
c23 = c.c23(index);
c33 = c.c33(index);

drr = dr.*dr;
drrr = drr.*dr;
dzz = dz.*dz;
dzzz = dzz.*dz;

drz   = dr.*dz;
drzz  = dr.*dzz;
drzzz = dr.*dzzz;
drrz  = drr.*dz;
drrzz = drr.*dzz;
drrzzz = drr.*dzzz;
drrrz  = drrr.*dz;
drrrzz = drrr.*dzz;

psi1 = c00       + c10.*dr    + c20.*drr    + c30.*drrr   + ...
       c01.*dz   + c11.*drz   + c21.*drrz   + c31.*drrrz  + ...
       c02.*dzz  + c12.*drzz  + c22.*drrzz  + c32.*drrrzz + ...
       c03.*dzzz + c13.*drzzz + c23.*drrzzz + c33.*drrr.*dzzz;

dsdr1 = c10       + 2*c20.*dr    + 3*c30.*drr   + ...
        c11.*dz   + 2*c21.*drz   + 3*c31.*drrz  + ...
        c12.*dzz  + 2*c22.*drzz  + 3*c32.*drrzz + ...
        c13.*dzzz + 2*c23.*drzzz + 3*c33.*drrzzz;
dsdr1 = dsdr1/g.dR;

dsdz1 =   c01      +   c11.*dr   +   c21.*drr   +   c31.*drrr   + ...
        2*c02.*dz  + 2*c12.*drz  + 2*c22.*drrz  + 2*c32.*drrrz  + ...
        3*c03.*dzz + 3*c13.*drzz + 3*c23.*drrzz + 3*c33.*drrrzz;
    
dsdz1 = dsdz1./g.dZ;

R1_inv = 1./R1;
br1 = -dsdz1.*R1_inv;
bz1 =  dsdr1.*R1_inv;

psiN = (psi1*g.ip_sign - g.ssimag)./(g.ssibry-g.ssimag);
fpol = polyval(g.fpol_coeffs,psiN);

% toroidal field
bt1 = g.bcentr*g.rzero.*R1_inv;
ic = find(psiN <= 1.0);
if ~isempty(ic)
%     fpol(ic) = dbsval(psiN,3,g.xk,g.mw,g.fpol_coeffs_spline);
    bt1(ic) = fpol(ic).*R1_inv(ic);
end

Bout.br = br1;
Bout.bz = bz1;
Bout.bphi = bt1;

if ~isempty(iroff) | ~isempty(izoff) %#ok<*OR2>
    if ~nowarn & ~isempty(iroff)
        warning(['Point off grid in R: R = ',num2str(R1(iroff)),'. [Rmin,Rmax] = [',num2str(g.r(2)),',',num2str(g.r(g.mw-1)),']'])
    end
    if ~nowarn & ~isempty(izoff)
        warning(['Point off grid in Z: Z = ',num2str(Z1(izoff)),'. [Zmin,Zmax] = [',num2str(g.z(2)),',',num2str(g.z(g.mh-1)),']'])
    end
    if toroidal_off_grid
        if ~nowarn
            warning(['Point(s) off grid --> returning toroidal field = 1 there'])
        end
        ierr = 0; Bout.br(mask == 1) = 0; Bout.bz(mask == 1) = 0; Bout.bphi(mask == 1) = 1;
        return;
    else
        ierr = 1; Bout = [];
        return;
    end    
end
ierr = 0;


