function g = flip_g(g)


method = 'flip_at_z0';

switch lower(method)
    case 'flip_at_z0'       
        z0 = g.z;
        psirz0 = g.psirz;
    otherwise
        error('This method is not implemented:%s',method)
end

g.psirz = fliplr(g.psirz);
g.bdry(2,:) = -g.bdry(2,:);
g.zmaxis = -g.zmaxis;
g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);


gfilename = strcat(g.filename,'_flip');
fprintf('writing gfile : %s\n',gfilename)
write_gfile(g,gfilename);
