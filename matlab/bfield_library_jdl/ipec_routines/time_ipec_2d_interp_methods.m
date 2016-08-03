clearvars;

run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\low\gpec\';
ipec = open_ipec_field(run_path);

R = 1.2;
Z = 0.1;
P = 10*pi/180;


ntry = 1;
tic;
for i = 1:ntry
    interp2(ipec.eq.r(:,1),ipec.eq.z(1,:),ipec.eq.b_r.',R,Z)
end
fprintf('This took %f\n',toc)


[rr,zz] = meshgrid(ipec.eq.r(:,1),ipec.eq.z(1,:));
for i = 1:ntry
    interp2(rr,zz,ipec.eq.b_r.',R,Z)
end
fprintf('This took %f\n',toc)


dr_grid = ipec.eq.r(2,1) - ipec.eq.r(1,1);
dz_grid = ipec.eq.z(1,2) - ipec.eq.z(1,1);
ir = floor((R - ipec.eq.r(1,1))/dr_grid) + 1;
iz = floor((Z - ipec.eq.z(1,1))/dz_grid) + 1;

tic;
for i = 1:ntry
%     ir = floor(interp1(ipec.eq.r(:,1),1:ipec.eq.nr,R));
%     iz = floor(interp1(ipec.eq.z(1,:),1:ipec.eq.nz,Z));
    dr_grid = ipec.eq.r(ir+1,1) - ipec.eq.r(ir,1);
    dz_grid = ipec.eq.z(1,iz+1) - ipec.eq.z(1,iz);
    Q = ipec.eq.b_r(ir:ir+1,iz:iz+1);
    drvec = [ipec.eq.r(ir+1,1) - R,R - ipec.eq.r(ir,1)];
    dzvec = [ipec.eq.z(1,iz+1) - Z;Z - ipec.eq.z(1,iz)];
    drvec*Q*dzvec/(dr_grid*dz_grid)
end
fprintf('This took %f\n',toc)
asdfdsafa
tic;
for i = 1:ntry
    ir = floor(interp1(ipec.eq.r(:,1),1:ipec.eq.nr,R));
    iz = floor(interp1(ipec.eq.z(1,:),1:ipec.eq.nz,Z));
    dr_grid = ipec.eq.r(ir+1,1) - ipec.eq.r(ir,1);
    dz_grid = ipec.eq.z(1,iz+1) - ipec.eq.z(1,iz);
    dr2 = ipec.eq.r(ir+1,1) - R;
    dr1 = R - ipec.eq.r(ir,1);
    dz2 = ipec.eq.z(1,iz+1) - Z;
    dz1 = Z - ipec.eq.z(1,iz);
    Q11 = ipec.eq.b_r(ir  ,iz);
    Q12 = ipec.eq.b_r(ir  ,iz+1);
    Q21 = ipec.eq.b_r(ir+1,iz);
    Q22 = ipec.eq.b_r(ir+1,iz+1);
    
    1/(dr_grid*dz_grid)*(Q11*dr2*dz2 + Q21*dr1*dz2 + Q12*dr2*dz1 + Q22*dr1*dz1)
end
fprintf('This took %f\n',toc)