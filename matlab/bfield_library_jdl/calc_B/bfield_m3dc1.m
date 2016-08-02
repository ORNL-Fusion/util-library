function [br,bz,bt,ierr] = bfield_m3dc1(b,r,phi,z,scale)
% fields are scaled to 1kA, and this shot is 4
if nargin < 5
    scale = 4;
end
ierr = 0;
nowarn = 1;
add_eq_fields = 1;
return_eq_fields_off_grid = 1;  % Do not set to 1 unless you add option to get eq fields from gfile!!!!
toroidal_field_off_grid = 0;

if any(isnan(r)) | any(isnan(z) | any(isnan(phi)) )
    error('r or z or phi is NaN')
end

% % good_inds = find( r >= b.r(1) & r <= b.r(end) & z >= b.z(1) & z <= b.z(end));

% check for points off grid in R
iroff = find(r < b.r(1) | r > b.r(end));


% check for points off grid in Z
izoff = find(z < b.z(1) | z > b.z(end));


% nr_b = length(b.r);
% nz_b = length(b.z);
% ir = NaN(nr_b,1);
% iz = NaN(nz_b,1);

% ir(good_inds) = interp1(b.r,1:b.nr,r(good_inds)).';
% iz(good_inds) = interp1(b.z,1:b.nz,z(good_inds)).';

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

if add_eq_fields == 1
    br_eq = b.br_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.br_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.br_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.br_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    bz_eq = b.bz_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bz_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.bz_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bz_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    bt_eq = b.bt_eq(ir+iz*b.nr).*(1-tr).*(1-tz) + b.bt_eq(ir+1+iz*b.nr).*tr.*(1-tz) + b.bt_eq(ir+(iz+1)*b.nr).*(1-tr).*tz + b.bt_eq(ir+1+(iz+1)*b.nr).*tr.*tz;
    
    br = br + br_eq;
    bz = bz + bz_eq;
    bt = bt + bt_eq;
end

if any(bt == 0)
%     figure; hold on; box on;
    aa=find(bt==0);
%     plot(r(aa),z(aa),'c.')   
%     error('this is bad')
    
    iroff = [iroff,aa];
    mask(iroff) = 1;
end

if ~isempty(iroff) | ~isempty(izoff)
    if toroidal_field_off_grid == 1
        if ~nowarn
            warning(['Point(s) off grid --> returning toroidal field = 1 there'])
        end
        ierr = 0; br(mask == 1) = 0; bz(mask == 1) = 0; bt(mask == 1) = 1;
        return;
    end    
    if return_eq_fields_off_grid == 1
        if ~nowarn
            warning(['Point(s) off grid  --> returning eq fields from gfile there'])
        end
        [Bout,ierr] = bfield_geq_bicub(b.g,r(mask == 1),z(mask == 1),nowarn);
        br(mask == 1) = Bout.br;
        bt(mask == 1) = Bout.bphi;
        bz(mask == 1) = Bout.bz;
        return;
    end
    if ~nowarn
        warning(['Point(s) off grid']) %: R = ',num2str(R1),'. [Rmin,Rmax] = [',num2str(g.r(1)),',',num2str(g.r(g.mw)),']'])
    end
    ierr = 1; br = []; bz = []; bt = [];
    return;
end


% % 
% % r,z,phi
% % br,bz,bt
