function plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty,quiet)
% plot_gfile(g,psi_min,psi_max,npsi,newfig,con_linesty,quiet)
%
% g can be g struction from readg_g3d, or gfile_name

plot_psi = 'psiN'; % psiN or psi
nRefine = 2;

if ischar(g)
    g = readg_g3d(g);
end
if nargin < 2
    npsi = 50;
end

if nargin < 5
    newfig = 1;
end
if nargin < 6
    con_linesty = '-';
end

rPlot = g.r;
zPlot = g.z;
switch lower(plot_psi)
    case lower('psiN')
        psiPlot = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);
        sepVal = 1;
        if nargin < 2
            psi_min = 0.6;
            psi_max = 1.1;
        end
        cLabel = '\psi_N';
    case lower('psi')
        psiPlot = g.psirz;        
        sepVal = g.ssibry;
        if nargin < 2
            psi_min = 0.6*(g.ssibry-g.ssimag) + g.ssimag;
            psi_max = 1.1*(g.ssibry-g.ssimag) + g.ssimag;        
        end
        cLabel = '\psi';
    otherwise
        error('Bad value for plot_psi')
end

%% Refine for better countours
for i = 1:nRefine
    [rPlot,zPlot,psiPlot] = refine_psi(rPlot,zPlot,2,g);
end



%%
if newfig == 1
    figure; hold on; box on;
end
% hh=contour(g.r,g.z,psiN_g.',linspace(psi_min,psi_max,npsi),con_linesty,'linewidth',lw);
% % set(hh,'linewidth',2)
% plot(lim_r,lim_z,'k.-')
% plot(bdry_r,bdry_z,'k--')

set(gcf,'color','w');
if isfield(g,'lim')
    plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),con_linesty,'linewidth',2)
end
contour(rPlot,zPlot,psiPlot.',linspace(psi_min,psi_max,npsi),'-','linewidth',1);
contour(rPlot,zPlot,psiPlot.',sepVal*[1,1],'k-','linewidth',2);
hc = colorbar;
set(hc,'fontsize',14);
colormap(turbo);
hc.Label.String = cLabel;
axis equal; axis tight;
ax = axis;
if isfield(g,'gfilename')
    h=text(ax(1)+(ax(2)-ax(1))*0.48,ax(3)+(ax(4)-ax(3))*0.97,g.gfilename,'fontsize',8);
    set(h,'interpreter','none');
end


% Find the xpt(s)
if isfield(g,'bdry')
    xpt_info = find_xpt_jl(g,1,1,1e-8,1);
    xr1 = xpt_info.rx;
    xz1 = xpt_info.zx;
    xr2 = xpt_info.rx2;
    xz2 = xpt_info.zx2;

    contour(rPlot,zPlot,psiPlot.',[1,1]*calc_psiN(g,xr2,xz2),'k-','linewidth',2)

    % contour(g.r,g.z,psiN_g.',[1,1]*2.25,'k-','linewidth',2)

    plot(xr1,xz1,'bx'); text(xr1+0.01,xz1,'x1','fontsize',8)
    plot(xr2,xz2,'b*'); text(xr2+0.02,xz2,'x2','fontsize',8)
end
% psi_x1 = g.ip_sign*get_psi_bicub(g,xr1,xz1);
% psi_x2 = g.ip_sign*get_psi_bicub(g,xr2,xz2);

% axis
rax = g.rmaxis;
zax = g.zmaxis;
plot(rax,zax,'bo')
xlabel('R [m]','fontsize',14)
ylabel('Z [m]','fontsize',14)
set(gca,'fontsize',14)
end

function [r2,z2,pn2] = refine_psi(r,z,fac,g)
r2 = linspace(r(2),r(end-1),length(r)*fac);
z2 = linspace(z(2),z(end-1),length(z)*fac);
for i = 1:length(z2)
    pn2(:,i) = calc_psiN(g,r2,z2(i)*ones(size(r2)));
end
end

