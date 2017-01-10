function [Br,Bz,Bphi]=bfield_grid(R,Z,P_rad,Bgrid,nowarn)


npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);


dr_grid = Bgrid.dR;
dz_grid = Bgrid.dZ;
dphi_grid = Bgrid.dphi;

for i = 1:npts
    
    Br(i) = 0;
    Bz(i) = 0;
    Bphi(i) = 0;
    
    P_rad_tmp = mod(P_rad(i)+2*pi/Bgrid.nsym,2*pi/Bgrid.nsym);
    
    if 0
        ir = floor(interp1(Bgrid.R(:),1:Bgrid.nr,R(i))); %stuck with interp as it gives NaN for OOB
        iz = floor(interp1(Bgrid.Z(:),1:Bgrid.nz,Z(i)));
        iphi = floor(interp1(Bgrid.phi,1:Bgrid.nphi,P_rad_tmp));
        if isnan(ir) || isnan(iz)
            if ~nowarn
                warning('Point(s) off grid in bfield_grid --> returning toroidal field = 1 there')
            end
            Bphi(i) = 1;
            continue;
        end
        if isnan(iphi)
            error('This should not happen because of mod')
        end
    else
        
        ir = floor( (R(i) - Bgrid.R(1))/dr_grid ) + 1;
        iz = floor( (Z(i) - Bgrid.Z(1))/dz_grid ) + 1;
        iphi = floor( (P_rad_tmp(i) - Bgrid.phi(1)) / dphi_grid ) + 1;
        if ir < 1 || ir >= Bgrid.nr || iz < 1 || iz >= Bgrid.nz
            warning('Point(s) off grid in bfield_grid --> returning toroidal field = 1 there')
            Bphi(i) = 1;
            continue;
        end
    end

    
    phi_fac = (P_rad_tmp - Bgrid.phi(iphi))/dphi_grid;
    
    dr2 = Bgrid.R(ir+1) - R(i);
    dr1 = dr_grid - dr2;
    dz2 = Bgrid.Z(iz+1) - Z(i);
    dz1 = dz_grid - dz2;
    
    QQ1 = Bgrid.Br(ir:ir+1,iz:iz+1,iphi);
    QQ2 = Bgrid.Br(ir:ir+1,iz:iz+1,iphi+1);
    Br(i) = ...
        (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
            phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        
    QQ1 = Bgrid.Bz(ir:ir+1,iz:iz+1,iphi);
    QQ2 = Bgrid.Bz(ir:ir+1,iz:iz+1,iphi+1);
    Bz(i) = ...
        (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
            phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        
    QQ1 = Bgrid.Bphi(ir:ir+1,iz:iz+1,iphi);
    QQ2 = Bgrid.Bphi(ir:ir+1,iz:iz+1,iphi+1);
    Bphi(i) = ...
        (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
            phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
end

% if nargout > 3
%     Btot = sqrt(Br.^2 + Bz.^2 + Bphi.^2);
% end

