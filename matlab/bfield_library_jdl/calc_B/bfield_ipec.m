function [Br,Bz,Bphi]=bfield_ipec(R,Z,P_rad,ipec,nowarn,field_choice)
% field_choice == 1: equilibrium only
% field_choice == 2: vacuum

npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);


for i = 1:npts
    
    ir = floor(interp1(ipec.eq.r(:,1),1:ipec.eq.nr,R(i)));
    iz = floor(interp1(ipec.eq.z(1,:),1:ipec.eq.nz,Z(i)));
    if isnan(ir) || isnan(iz)
        if ~nowarn
            warning(['Point(s) off grid in bfield_ipec --> returning toroidal field = 1 there'])
        end
        Br = 0;
        Bz = 0;
        Bphi = 1;
        return;
    end
    dr_grid = ipec.eq.r(ir+1,1) - ipec.eq.r(ir,1);
    dz_grid = ipec.eq.z(1,iz+1) - ipec.eq.z(1,iz);

    dr2 = ipec.eq.r(ir+1,1) - R(i);
    dr1 = dr_grid - dr2;
    dz2 = ipec.eq.z(1,iz+1) - Z(i);
    dz1 = dz_grid - dz2;

    QQ = ipec.eq.b_r(ir:ir+1,iz:iz+1);
    Br(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
    QQ = ipec.eq.b_z(ir:ir+1,iz:iz+1);
    Bz(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
    QQ = ipec.eq.b_phi(ir:ir+1,iz:iz+1);
    Bphi(i) = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);

    if field_choice == 2
        QQ = ipec.vac.rb_r(ir:ir+1,iz:iz+1);
        rbr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ib_r(ir:ir+1,iz:iz+1);
        ibr = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Br(i) = Br(i) + rbr*cos(ipec.vac.n*P_rad) + ibr*sin(ipec.vac.n*P_rad);

        QQ = ipec.vac.rb_z(ir:ir+1,iz:iz+1);
        rbz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ib_z(ir:ir+1,iz:iz+1);
        ibz = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Bz(i) = Bz(i) + rbz*cos(ipec.vac.n*P_rad) + ibz*sin(ipec.vac.n*P_rad);      
        
        QQ = ipec.vac.rb_phi(ir:ir+1,iz:iz+1);
        rbphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
        QQ = ipec.vac.ib_phi(ir:ir+1,iz:iz+1);
        ibphi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);        
        Bphi(i) = Bphi(i) + rbphi*cos(ipec.vac.n*P_rad) + ibphi*sin(ipec.vac.n*P_rad);      
    end
        
end