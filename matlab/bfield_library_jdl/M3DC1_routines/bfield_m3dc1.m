function [br,bz,bt,ierr] = bfield_m3dc1(b,r,phi,z,scale)
% fields are scaled to 1kA, and this shot is 4
if nargin < 5
    scale = 4;
end
ierr = 0;
nowarn = 0;
add_eq_fields = 1;
return_eq_fields_off_grid = 0;

if any(isnan(r)) | any(isnan(z) | any(isnan(phi)) )
    error('r or z or phi is NaN')
end

% check for points off grid 
iroff = find(r < b.r(1) | r > b.r(end));
if ~isempty(iroff)
    if ~nowarn
        if return_eq_fields_off_grid == 1
            warning(['Point(s) off grid in R --> returning eq fields there'])
        else
            warning(['Point(s) off grid in R']) %: R = ',num2str(R1),'. [Rmin,Rmax] = [',num2str(g.r(1)),',',num2str(g.r(g.mw)),']'])
        end
    end
    if ~return_eq_fields_off_grid
        ierr = 1; br = []; bz = []; bt = [];
        return;
    end
end
izoff = find(z < b.z(1) | z > b.z(end));
if ~isempty(izoff)
    if ~nowarn
        if return_eq_fields_off_grid == 1
            warning(['Point(s) off grid in Z --> returning eq fields there'])
        else
            warning(['Point(s) off grid in Z']) %: Z = ',num2str(Z1),'. [Zmin,Zmax] = [',num2str(g.z(1)),',',num2str(g.z(g.mh)),']'])
        end
    end
    if ~return_eq_fields_off_grid
        ierr = 1; br = []; bz = []; bt = [];
        return;
    end
end
ir = interp1(b.r,1:b.nr,r).';
iz = interp1(b.z,1:b.nz,z).';

mask = zeros(size(r));
mask(iroff) = 1;
mask(izoff) = 1;
ir(mask==1) = 1;
iz(mask==1) = 1;

tr = ir - floor(ir);
tz = iz - floor(iz);
ir = floor(ir);
iz = floor(iz);

br_amp = b.br_amp(ir+iz*b.nr).*(1-tr).*(1-tz) + b.br_amp(ir+1+iz*b.nr).*tr.*(1-tz) + b.br_amp(ir+(iz+1)*b.nr).*(1-tr).*tz + b.br_amp(ir+1+(iz+1)*b.nr).*tr.*tz;
bz_amp = b.bz_amp(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bz_amp(ir+1+iz*b.nr).*tr.*(1-tz) + b.bz_amp(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bz_amp(ir+1+(iz+1)*b.nr).*tr.*tz;
bt_amp = b.bt_amp(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bt_amp(ir+1+iz*b.nr).*tr.*(1-tz) + b.bt_amp(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bt_amp(ir+1+(iz+1)*b.nr).*tr.*tz;
br_phase = b.br_phase(ir+iz*b.nr).*(1-tr).*(1-tz) + b.br_phase(ir+1+iz*b.nr).*tr.*(1-tz) + b.br_phase(ir+(iz+1).*b.nr).*(1-tr).*tz + b.br_phase(ir+1+(iz+1)*b.nr).*tr.*tz;
bz_phase = b.bz_phase(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bz_phase(ir+1+iz*b.nr).*tr.*(1-tz) + b.bz_phase(ir+(iz+1).*b.nr).*(1-tr).*tz + b.bz_phase(ir+1+(iz+1)*b.nr).*tr.*tz;
bt_phase = b.bt_phase(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bt_phase(ir+1+iz*b.nr).*tr.*(1-tz) + b.bt_phase(ir+(iz+1).*b.nr).*(1-tr).*tz + b.bt_phase(ir+1+(iz+1)*b.nr).*tr.*tz;

br = scale*br_amp.*sin(br_phase+phi.'*b.nn);
bz = scale*bz_amp.*sin(bz_phase+phi.'*b.nn);
bt = scale*bt_amp.*sin(bt_phase+phi.'*b.nn);

br(mask == 1) = 0;
bt(mask == 1) = 0;
bz(mask == 1) = 0;

if add_eq_fields
    br_eq = b.br_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.br_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.br_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.br_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    bz_eq = b.bz_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bz_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.bz_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bz_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    bt_eq = b.bt_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bt_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.bt_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bt_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    
    br = br + br_eq;
    bz = bz + bz_eq;
    bt = bt + bt_eq;
end

