function plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty)
% plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty)

% gfile_name = 'C:\Work\EMC3_revival\gfiles\g140508.00403';
% gfile_name = 'C:\Work\EMC3_revival\gfiles\g1090814013.00743_604';
% g = readg_g3d(gfile_name);

if nargin < 2
    npsi = 50;
    psi_min = 0.6;
    psi_max = 1.1;
end

if nargin < 5
    newfig = 1;
end
if nargin < 6
    con_linesty = '-';
end

lw = 2;

lim_r = g.lim(1,g.lim(1,:) > 0);
lim_z = g.lim(2,g.lim(1,:) > 0);
psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);
bdry_r = g.bdry(1,:);
bdry_z = g.bdry(2,:);

if newfig == 1
    figure; hold on; box on;
end
hh=contour(g.r,g.z,psiN_g.',linspace(psi_min,psi_max,npsi),con_linesty,'linewidth',lw);
% set(hh,'linewidth',2)
plot(lim_r,lim_z,'k.-')
plot(bdry_r,bdry_z,'k--')

% Find the xpt(s)
% xpt_info = find_xpt_jl(g,1,1,1e-8,1,[],[],0.1);
xpt_info = find_xpt_jl(g,1,1,1e-8,0);
xr1 = xpt_info.rx;
xz1 = xpt_info.zx;
xr2 = xpt_info.rx2;
xz2 = xpt_info.zx2;

plot(xr1,xz1,'bx')
plot(xr2,xz2,'bx')

psi_x1 = g.ip_sign*get_psi_bicub(g,xr1,xz1);
psi_x2 = g.ip_sign*get_psi_bicub(g,xr2,xz2);

% axis
rax = g.rmaxis;
zax = g.zmaxis;
plot(rax,zax,'bx')
xlabel('R [m]','fontsize',14)
ylabel('Z [m]','fontsize',14)
set(gca,'fontsize',14)

