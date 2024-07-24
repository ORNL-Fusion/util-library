function g = symm_g(g,method)

if nargin < 2
    method = 'avg_zax';
end

if mod(g.mh,2) == 0
    error('Assumed that g.mh is odd')
end

switch lower(method)
    case 'avg_z0'       
        z0 = g.z;
        psirz0 = g.psirz;
    case 'avg_zax'
        z0min = min(abs([g.z(end) + g.zmaxis,g.z(end),g.z(1),g.z(1)+g.zmaxis]));
        z0 = linspace(-z0min,z0min,g.mh);
        psirz0 = zeros(g.mw,g.mh);
        for i = 1:g.mw
            psirz0(i,:) = calc_psi(g,g.r(i)*ones(size(z0)),z0 + g.zmaxis)./g.ip_sign;
        end
end

psirz = zeros(g.mw, g.mh);

% Loop through the first half of Z to symmetrize the matrix
for j = 1:floor(g.mh / 2)
    for i = 1:g.mw
        psirz(i, j) = 0.5 * (psirz0(i, j) + psirz0(i, g.mh - j + 1));
        psirz(i, g.mh - j + 1) = psirz(i, j);
    end
end

jmid = ceil(g.mh / 2);
for i = 1:g.mw
%     psirz(i, jmid) = 0.5 * (psirz0(i, jmid - 1) + psirz0(i, jmid + 1));
    psirz(i, jmid) = psirz0(i, jmid);
end


g.psirz = psirz;
g.z = z0;
g.zmaxis = 0;

g.zdim = g.z(end) - g.z(1);

g.dZ = g.zdim/(g.mh-1);

g.z = zeros(1,g.mh);

for i = 0:g.mh-1
    g.z(1,i+1) = g.zmid - 0.5*g.zdim + g.dZ*i;
end

g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);


gfilename = strcat(g.filename,'_symm');
fprintf('writing gfile : %s\n',gfilename)
write_gfile(g,gfilename);
