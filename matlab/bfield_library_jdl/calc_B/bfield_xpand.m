function [Br,Bz,Bphi,Btot]=bfield_xpand(R,Z,P_rad,field,nowarn,field_choice)
% field_choice == 0: total
% field_choice == 1: vacuum

npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);


for i = 1:npts
    
    Br(i) = 0;
    Bz(i) = 0;
    Bphi(i) = 0;
    
    ir = floor(interp1(field.r(:),1:field.nr,R(i))); %stuck with interp as it gives NaN for OOB
    iz = floor(interp1(field.z(:),1:field.nz,Z(i)));
    P_rad(i) = mod(P_rad(i)+2*pi,2*pi);
    iphi = floor(interp1(field.phi,1:field.nphi,P_rad(i)));
    error('Replace above because interp1 is hella slow: check test_bfield_grid')
%     iphi = floor( (P_rad(i) - field.phi(1)) / (field.phi(2)-field.phi(1)) ) + 1;
    
    if isnan(ir) || isnan(iz)
        if ~nowarn
            warning(['Point(s) off grid in bfield_xpand --> returning toroidal field = 1 there'])
        end
        Bphi(i) = 1;
        continue;
    end
    if isnan(iphi)
        error('This should not happen because of mod')
    end
    
    dr_grid = field.r(ir+1) - field.r(ir);
    dz_grid = field.z(iz+1) - field.z(iz);
    dphi_grid = field.phi(iphi+1) - field.phi(iphi);   
    phi_fac = (P_rad(i) - field.phi(iphi))/dphi_grid;
       
    dr2 = field.r(ir+1) - R(i);
    dr1 = dr_grid - dr2;
    dz2 = field.z(iz+1) - Z(i);
    dz1 = dz_grid - dz2;    

    if field_choice == 0
        QQ1 = field.Br(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Br(ir:ir+1,iz:iz+1,iphi+1);        
        Br(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ1 = field.Bz(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Bz(ir:ir+1,iz:iz+1,iphi+1);        
        Bz(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        QQ1 = field.Bphi(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Bphi(ir:ir+1,iz:iz+1,iphi+1);        
        Bphi(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
    elseif field_choice == 1
        QQ1 = field.Brvac(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Brvac(ir:ir+1,iz:iz+1,iphi+1);        
        Br(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ1 = field.Bzvac(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Bzvac(ir:ir+1,iz:iz+1,iphi+1);        
        Bz(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        QQ1 = field.Bphivac(ir:ir+1,iz:iz+1,iphi);
        QQ2 = field.Bphivac(ir:ir+1,iz:iz+1,iphi+1);        
        Bphi(i) = ...
            (1-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1)/(dr_grid*dz_grid) + ...
                phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        
%         % Need to add EFIT?? 
%         Bout = bfield_geq_bicub(field.g,R(i),Z(i),nowarn);
%         Br(i) = Br(i) + Bout.br;
%         Bphi(i) = Bphi(i) + Bout.bphi;
%         Bz(i) = Bz(i) + Bout.bz;
    else
        error('bad field_choice')
    end
       
end

if nargout > 3
    Btot = sqrt(Br.^2 + Bz.^2 + Bphi.^2);
end