function [Br,Bz,Bphi,Btot]=bfield_ipec(R,Z,P_rad,ipec,nowarn,field_choice)
% field_choice == 0: equilibrium only
% field_choice == 1: vacuum + eq
% field_choice == 2: pert + eq
% field_choice == 3: vacuum pert only
% field_choice == 4: pert only

npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);


for i = 1:npts
    
    Br(i) = 0;
    Bz(i) = 0;
    Bphi(i) = 0;
    
    ir = floor(interp1(ipec.eq.r(:,1),1:ipec.eq.nr,R(i))); %stuck with interp as it gives NaN for OOB
    iz = floor(interp1(ipec.eq.z(1,:),1:ipec.eq.nz,Z(i)));
    error('Replace above because interp1 is hella slow: check test_bfield_grid')
    if isnan(ir) || isnan(iz)
        if ~nowarn
            warning(['Point(s) off grid in bfield_ipec --> returning toroidal field = 1 there'])
        end
        Bphi(i) = 1;
        continue;
    end
    dr_grid = ipec.eq.r(ir+1,1) - ipec.eq.r(ir,1);
    dz_grid = ipec.eq.z(1,iz+1) - ipec.eq.z(1,iz);

    dr2 = ipec.eq.r(ir+1,1) - R(i);
    dr1 = dr_grid - dr2;
    dz2 = ipec.eq.z(1,iz+1) - Z(i);
    dz1 = dz_grid - dz2;

    if field_choice < 3
        QQ = ipec.eq.br(ir:ir+1,iz:iz+1);
        Br(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.eq.bz(ir:ir+1,iz:iz+1);
        Bz(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.eq.bphi(ir:ir+1,iz:iz+1);
        Bphi(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
    end
    if field_choice == 0
        continue;
    elseif (field_choice == 1 || field_choice == 3)
        sinphi = sin(ipec.vac.n*P_rad(i));
        cosphi = cos(ipec.vac.n*P_rad(i));
        
        QQ = ipec.vac.rbr(ir:ir+1,iz:iz+1);
        rbr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ibr(ir:ir+1,iz:iz+1);
        ibr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Br(i) = Br(i) + rbr*cosphi + ibr*sinphi;

        QQ = ipec.vac.rbz(ir:ir+1,iz:iz+1);
        rbz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ibz(ir:ir+1,iz:iz+1);
        ibz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Bz(i) = Bz(i) + rbz*cosphi + ibz*sinphi;      
        
        QQ = ipec.vac.rbphi(ir:ir+1,iz:iz+1);
        rbphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ibphi(ir:ir+1,iz:iz+1);
        ibphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        Bphi(i) = Bphi(i) + rbphi*cosphi + ibphi*sinphi;
    elseif (field_choice == 2 || field_choice == 4)
        sinphi = sin(ipec.pert.n*P_rad(i));
        cosphi = cos(ipec.pert.n*P_rad(i));
        
        QQ = ipec.pert.rbr(ir:ir+1,iz:iz+1);
        rbr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.pert.ibr(ir:ir+1,iz:iz+1);
        ibr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        Br(i) = Br(i) + rbr*cosphi + ibr*sinphi;

        QQ = ipec.pert.rbz(ir:ir+1,iz:iz+1);
        rbz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.pert.ibz(ir:ir+1,iz:iz+1);
        ibz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Bz(i) = Bz(i) + rbz*cosphi + ibz*sinphi;      
        
        QQ = ipec.pert.rbphi(ir:ir+1,iz:iz+1);
        rbphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.pert.ibphi(ir:ir+1,iz:iz+1);
        ibphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Bphi(i) = Bphi(i) + rbphi*cosphi + ibphi*sinphi;          
    else
        error('Bad field choice')
    end
end

if nargout > 3
    Btot = sqrt(Br.^2 + Bz.^2 + Bphi.^2);
end