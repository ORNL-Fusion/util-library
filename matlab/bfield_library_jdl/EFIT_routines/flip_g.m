function g = flip_g(g)
% Routine to flip the gfile flux 
% Only flip about Z=0 implemented so far
% JD Lore 2025

method = 'flip_at_z0';

switch lower(method)
    case 'flip_at_z0'       
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
