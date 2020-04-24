function plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty,quiet)
% plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty)

% gfile_name = 'C:\Work\EMC3_revival\gfiles\g140508.00403';
% gfile_name = 'C:\Work\EMC3_revival\gfiles\g1090814013.00743_604';
% g = readg_g3d(gfile_name);

if isstr(g)
    g = readg_g3d(g);
end
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

% lw = 2;
% 
% lim_r = g.lim(1,g.lim(1,:) > 0);
% lim_z = g.lim(2,g.lim(1,:) > 0);
psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);
% bdry_r = g.bdry(1,:);
% bdry_z = g.bdry(2,:);

if newfig == 1
    figure; hold on; box on;
end
% hh=contour(g.r,g.z,psiN_g.',linspace(psi_min,psi_max,npsi),con_linesty,'linewidth',lw);
% % set(hh,'linewidth',2)
% plot(lim_r,lim_z,'k.-')
% plot(bdry_r,bdry_z,'k--')


% figure; hold on; box on;
if isfield(g,'lim')
    plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),con_linesty,'linewidth',2)
end
contour(g.r,g.z,psiN_g.',linspace(psi_min,psi_max,npsi),'-','linewidth',1);
contour(g.r,g.z,psiN_g.',[1,1],'k-','linewidth',2);
hc = colorbar;
set(hc,'fontsize',14);
colormap(colorflipper(512,'jet'));
axis equal; axis tight;
% axis([1,2.4,-1.5,1.5])
ax = axis;
if isfield(g,'gfilename')
    h=text(ax(1)+(ax(2)-ax(1))*0.48,ax(3)+(ax(4)-ax(3))*0.97,g.gfilename,'fontsize',8);
    set(h,'interpreter','none');
end



% Find the xpt(s)
xpt_info = find_xpt_jl(g,1,1,1e-8,1);
xr1 = xpt_info.rx;
xz1 = xpt_info.zx;
xr2 = xpt_info.rx2;
xz2 = xpt_info.zx2;

contour(g.r,g.z,psiN_g.',[1,1]*calc_psiN(g,xr2,xz2),'k-','linewidth',2)

% contour(g.r,g.z,psiN_g.',[1,1]*2.25,'k-','linewidth',2)

plot(xr1,xz1,'bx'); text(xr1+0.01,xz1,'x1','fontsize',8)
plot(xr2,xz2,'b*'); text(xr2+0.02,xz2,'x2','fontsize',8)

% psi_x1 = g.ip_sign*get_psi_bicub(g,xr1,xz1);
% psi_x2 = g.ip_sign*get_psi_bicub(g,xr2,xz2);

% axis
rax = g.rmaxis;
zax = g.zmaxis;
plot(rax,zax,'bo')
xlabel('R [m]','fontsize',14)
ylabel('Z [m]','fontsize',14)
set(gca,'fontsize',14)

