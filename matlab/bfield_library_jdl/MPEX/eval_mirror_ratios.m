function ratios = eval_mirror_ratios(bfield,geo)
if nargin < 2
    [geo] = get_Proto_geometry;
end


z_3p5 = (geo.cmin(4) - geo.cmax(3)) + geo.cmax(3);
r_3p5 = 1e-3;

z_1 = (geo.cmin(1) + geo.cmax(1))/2;
r_1 = 1e-3;

z_2 = (geo.cmin(2) + geo.cmax(2))/2;
r_2 = 1e-3;

z_5 = (geo.cmin(5) + geo.cmax(5))/2;
r_5 = 1e-3;

z_6 = (geo.cmin(6) + geo.cmax(6))/2;
r_6 = 1e-3;

phieval = 0;
switch bfield.type   
    case 'just_coils'        
        [~,~,~,B_3p5]=bfield_bs_cyl(r_3p5,phieval,z_3p5,bfield.coil,bfield.current);
        [~,~,~,B_1]=bfield_bs_cyl(r_1,phieval,z_1,bfield.coil,bfield.current);
        [~,~,~,B_2]=bfield_bs_cyl(r_2,phieval,z_2,bfield.coil,bfield.current);
        [~,~,~,B_5]=bfield_bs_cyl(r_5,phieval,z_5,bfield.coil,bfield.current);
        [~,~,~,B_6]=bfield_bs_cyl(r_6,phieval,z_6,bfield.coil,bfield.current);        
    otherwise
        error('did not recognize bfield.type: ',bfield.type)
end

ratios.R_u1 = B_1/B_3p5;
ratios.R_u2 = B_2/B_3p5;
ratios.R_d5 = B_5/B_3p5;
ratios.R_d6 = B_6/B_3p5;