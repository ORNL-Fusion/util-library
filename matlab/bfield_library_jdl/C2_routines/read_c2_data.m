function data = read_c2_data(run_path,ndomain,time_str)
% time_str indicates the "000001" of files like 2d.000001.003.dat,
% and can also be "final".

for idomain = 1:ndomain
    istr_app = '';
    if idomain < 100 && idomain > 9
        istr_app = '0';
    elseif idomain < 10
        istr_app = '00';
    end
    
    data_fname = [run_path,'2d.',time_str,'.',istr_app,num2str(idomain),'.dat'];
    
    fid = fopen(data_fname);
    node_num = fscanf(fid,'%d',2);
    num_pts = prod(node_num);
    nx = node_num(1);
    ny = node_num(2);
    
    data_tmp = fscanf(fid,'%d %d %e %e %e %e %e %e %e %e',10*num_pts);
    data_tmp = reshape(data_tmp,10,num_pts);

    xc1d = data_tmp(3,:);
    yc1d = data_tmp(4,:);
    Te1d = data_tmp(5,:);
    Ti1d = data_tmp(6,:);
    np1d = data_tmp(7,:);
    nn1d = data_tmp(8,:);
    u1d  = data_tmp(9,:);
    bx1d = data_tmp(10,:);
    
    data.xc2d{idomain} = reshape(xc1d,ny,nx);
    data.yc2d{idomain} = reshape(yc1d,ny,nx);
    data.Te2d{idomain} = reshape(Te1d,ny,nx);
    data.Ti2d{idomain} = reshape(Ti1d,ny,nx);
    data.np2d{idomain} = reshape(np1d,ny,nx);
    data.nn2d{idomain} = reshape(nn1d,ny,nx);
    data.u2d{idomain}  = reshape(u1d,ny,nx);
    data.bx2d{idomain} = reshape(bx1d,ny,nx);
    
    fclose(fid);
end
