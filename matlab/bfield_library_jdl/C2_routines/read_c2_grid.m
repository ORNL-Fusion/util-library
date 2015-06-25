function grid = read_c2_grid(grid_dir,ndomain)

for idomain = 1:ndomain
    istr_app = '';
    if idomain < 100 && idomain > 9
        istr_app = '0';
    elseif idomain < 10
        istr_app = '00';
    end
    
    grid_fname = [grid_dir,'grid.',istr_app,num2str(idomain),'.dat'];
    
    fid = fopen(grid_fname);
    data = fscanf(fid,'%s',1);
    if isempty(strfind((lower(data)),'domain')); error('output not in assumed order (domain)'); end
    domain_id = fscanf(fid,'%d',1);
    
    data = fscanf(fid,'%s',1);
    if isempty(strfind((lower(data)),'nodenum')); error('output not in assumed order (nodenum)'); end
    node_num = fscanf(fid,'%d',2);
    num_pts = prod(node_num);
    nx = node_num(1);
    ny = node_num(2);
    
    data = fscanf(fid,'%s',1);
    if isempty(strfind((lower(data)),'connect')); error('output not in assumed order (connect)'); end
    connect = fscanf(fid,'%d',4);
    
    data = fscanf(fid,'%s',1);
    if isempty(strfind((lower(data)),'data')); error('output not in assumed order (data)'); end
    grid_data = fscanf(fid,'%e',num_pts*5);
    grid_data = reshape(grid_data,5,num_pts);
    x1d = grid_data(1,:);
    y1d = grid_data(2,:);
    x2d = reshape(x1d,ny,nx);
    y2d = reshape(y1d,ny,nx);
    
    grid.x2d{idomain} = x2d;
    grid.y2d{idomain} = y2d;
    fclose(fid);
    
    grid.nx(idomain) = nx;
    grid.ny(idomain) = ny;
    grid.num_pts(idomain) = num_pts;
    
end
grid.ndomain = ndomain;

grid.Zoffsets(1) = 0;
for iz = 2:grid.ndomain
    grid.Zoffsets(iz) = grid.Zoffsets(iz-1) + (grid.ny(iz-1) - 1)*(grid.nx(iz-1) - 1);
end
% iz = grid.nz+1;
% Zoff_all = grid.Zoffsets(iz-1) + (grid.nr_array(iz-1) - 1)*(grid.np_array(iz-1) - 1);

