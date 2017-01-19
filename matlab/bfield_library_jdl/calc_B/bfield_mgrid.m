function [Br,Bz,Bphi]=bfield_mgrid(R,Z,P_rad,mgrid,nowarn)

error('I just wrote the coil loop into bfield_grid')

nc = length(mgrid.raw_coil_cur);

Br_c = zeros(nc,1);
Bz_c = zeros(nc,1);
Bphi_c = zeros(nc,1);
for i = 1:nc
    scale_fac = mgrid.scale_factor(i);
    if abs(scale_fac) > 1e-10
        [Br_c(i),Bz_c(i),Bphi_c(i)]=bfield_grid(R,Z,P_rad,mgrid,nowarn);
    else
        Br_c(i) = 0;
        Bz_c(i) = 0;
        Bphi_c(i) = 0;
    end
end

Br = sum(Br_c);
Bz = sum(Bz_c);
Bphi = sum(Bphi_c);