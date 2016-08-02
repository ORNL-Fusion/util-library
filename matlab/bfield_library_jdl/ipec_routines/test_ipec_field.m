clearvars;

run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\low\gpec\';
ipec = open_ipec_field(run_path);

R = 1.2;
Z = 0.1;
P = 10*pi/180;



ir = floor(interp1(ipec.eq.r(:,1),1:ipec.eq.nr,R));
iz = floor(interp1(ipec.eq.z(1,:),1:ipec.eq.nz,Z));
dr_grid = ipec.eq.r(ir+1,1) - ipec.eq.r(ir,1);
dz_grid = ipec.eq.z(1,iz+1) - ipec.eq.z(1,iz);
dr2 = ipec.eq.r(ir+1,1) - R;
dr1 = dr_grid - dr2;
dz2 = ipec.eq.z(1,iz+1) - Z;
dz1 = dz_grid - dz2;


QQ = ipec.eq.b_r(ir:ir+1,iz:iz+1);
b_r = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
QQ = ipec.eq.b_z(ir:ir+1,iz:iz+1);
b_z = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);
QQ = ipec.eq.b_phi(ir:ir+1,iz:iz+1);
b_phi = (QQ(1,1)*dr2*dz2 + QQ(2,1)*dr1*dz2 + QQ(1,2)*dr2*dz1 + QQ(2,2)*dr1*dz1)/(dr_grid*dz_grid);