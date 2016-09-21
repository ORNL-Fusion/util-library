function [Bout,ierr] = bfield_geq_bicub(g,R1,Z1,nowarn)
%[Bout,ierr] = bfield_geq_bicub(g,R1,Z1,nowarn)
if nargin < 4
    nowarn = 0;
end

if iscolumn(R1)
    R1 = R1.';
    Z1 = Z1.';
end

toroidal_off_grid = 0;
if isfield(g,'toroidal_off_grid')
    toroidal_off_grid = g.toroidal_off_grid;
end        

ir = floor( (R1 - g.r(1))./g.dR ) + 1;
iz = floor( (Z1 - g.z(1))./g.dZ ) + 1;

if any(isnan(R1)) | any(isnan(Z1))
    error('R1 or Z1 is NaN')
end

% check for points off grid
iroff = find(ir <= 2 | ir >= g.mw - 1);
izoff = find(iz <= 1 | iz >= g.mh - 1);

mask = zeros(size(R1));
mask(iroff) = 1;
mask(izoff) = 1;
ir(mask==1) = 3;  % just a value that will not error
iz(mask==1) = 3;

dir(:,1) = (R1 - g.r(ir))./g.dR;
diz(:,1) = (Z1 - g.z(iz))./g.dZ;

dir2 = dir.*dir;
dir3 = dir2.*dir;
diz2 = diz.*diz;
diz3 = diz2.*diz;

index = iz + g.mh*(ir-1);
c_bi = g.bicub_coeffs(index,:,:);

psi1 = c_bi(:,1,1)        + c_bi(:,2,1).*dir        + c_bi(:,3,1).*dir2        + c_bi(:,4,1).*dir3       + ...
       c_bi(:,1,2).*diz   + c_bi(:,2,2).*dir.*diz   + c_bi(:,3,2).*dir2.*diz   + c_bi(:,4,2).*dir3.*diz  + ...
       c_bi(:,1,3).*diz2  + c_bi(:,2,3).*dir.*diz2  + c_bi(:,3,3).*dir2.*diz2  + c_bi(:,4,3).*dir3.*diz2 + ...
       c_bi(:,1,4).*diz3  + c_bi(:,2,4).*dir.*diz3  + c_bi(:,3,4).*dir2.*diz3  + c_bi(:,4,4).*dir3.*diz3;

dsdr1 = c_bi(:,2,1)        + 2*c_bi(:,3,1).*dir        + 3*c_bi(:,4,1).*dir2       + ...
        c_bi(:,2,2).*diz   + 2*c_bi(:,3,2).*dir.*diz   + 3*c_bi(:,4,2).*dir2.*diz  + ...
        c_bi(:,2,3).*diz2  + 2*c_bi(:,3,3).*dir.*diz2  + 3*c_bi(:,4,3).*dir2.*diz2 + ...
        c_bi(:,2,4).*diz3  + 2*c_bi(:,3,4).*dir.*diz3  + 3*c_bi(:,4,4).*dir2.*diz3;
dsdr1 = dsdr1/g.dR;

dsdz1 =   c_bi(:,1,2)       +   c_bi(:,2,2).*dir       +   c_bi(:,3,2).*dir2       +   c_bi(:,4,2).*dir3       + ...
        2*c_bi(:,1,3).*diz  + 2*c_bi(:,2,3).*dir.*diz  + 2*c_bi(:,3,3).*dir2.*diz  + 2*c_bi(:,4,3).*dir3.*diz  + ...
        3*c_bi(:,1,4).*diz2 + 3*c_bi(:,2,4).*dir.*diz2 + 3*c_bi(:,3,4).*dir2.*diz2 + 3*c_bi(:,4,4).*dir3.*diz2;
dsdz1 = dsdz1./g.dZ;

br1 = -dsdz1./R1.';
bz1 =  dsdr1./R1.';

psiN = (psi1*g.ip_sign - g.ssimag)./(g.ssibry-g.ssimag);
fpol = polyval(g.fpol_coeffs,psiN);

% toroidal field
bt1 = g.bcentr*g.rzero./R1;
ic = find(psiN <= 1.0);
if ~isempty(ic)
    bt1(ic) = fpol(ic)./R1(ic).';
end

Bout.br = br1;
Bout.bz = bz1;
Bout.bphi = bt1.';

if ~isempty(iroff) | ~isempty(izoff) %#ok<*OR2>
    if ~nowarn & ~isempty(iroff)
        warning(['Point off grid in R: R = ',num2str(R1(iroff)),'. [Rmin,Rmax] = [',num2str(g.r(1)),',',num2str(g.r(g.mw)),']'])
    end
    if ~nowarn & ~isempty(izoff)
        warning(['Point off grid in Z: Z = ',num2str(Z1(izoff)),'. [Zmin,Zmax] = [',num2str(g.z(1)),',',num2str(g.z(g.mh)),']'])
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


