function plot_gfile(g,varargin)
% plot_gfile(g)
% plot_gfile(g,'Name',value,...)
%
% plot_psi is 'psi' (default) or 'psiN' 
% plot_b_components plots Br, Bpol, Btor when set to 1
%
% g can be g structure from readg_g3d, or gfile_name

nRefine = 1;

%% Handle inputs
psi_min = [];
psi_max = [];
npsi = 50;
newfig = 1;
line_style = '-';
quiet = 1;
debug_bpol_est = 0;
plot_psi = 'psi';
plot_b_components = 0;

if mod(length(varargin),2) ~= 0
    error('Optional inputs must be name/value pairs')
end

for i = 1:2:length(varargin)
    if ~(ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i})))
        error('Optional input names must be character vectors or strings')
    end
    switch lower(char(varargin{i}))
        case {'psi_min','psimin'}
            psi_min = varargin{i+1};
        case {'psi_max','psimax'}
            psi_max = varargin{i+1};
        case 'npsi'
            npsi = varargin{i+1};
        case {'newfig','newfigure'}
            newfig = varargin{i+1};
        case {'linestyle','line_style'}
            line_style = varargin{i+1};
        case {'plot_psi','plotpsi'}
            plot_psi = varargin{i+1};
        case {'plot_b_components','plotbcomponents','plot_b'}
            plot_b_components = varargin{i+1};
        otherwise
            error('Unknown option: %s',char(varargin{i}))
    end
end

if ischar(g) || (isstring(g) && isscalar(g))
    g = readg_g3d(char(g),0,quiet);
end

%% Handle plot_psi switch
switch lower(plot_psi)
    case lower('psiN')
        psiPlot = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag); 
        sepVal = 1;
        if isempty(psi_min)
            psi_min = 0.6;
        end
        if isempty(psi_max)
            psi_max = 1.1;
        end
        
        cLabel = '\psi_N';
    case lower('psi')
        psiPlot = g.ip_sign*g.psirz; 
        sepVal = g.ip_sign*g.ssibry;
        axisVal = g.ip_sign*g.ssimag;
        if isempty(psi_min)
            psi_min = 0.6*(sepVal-axisVal) + axisVal;
        end
        if isempty(psi_max)
            psi_max = 1.1*(sepVal-axisVal) + axisVal;
        end
        cLabel = '\psi';
    otherwise
        error('Bad value for plot_psi')
end

%% Refine for better contours
if nRefine > 0
    % rPlot = g.r(1:end); zPlot = g.z(1:end);
    % rPlot = g.r(2:end-1); zPlot = g.z(2:end-1);
    [rPlot,zPlot,psiPlot] = refine_psi(2^nRefine,g,plot_psi);
else
    rPlot = g.r; zPlot = g.z;
    psiPlot = psiPlot.';
end
psiLevels = linspace(psi_min,psi_max,npsi);

%%
if newfig == 1
    figure; box on;
end
hold on;
set(gcf,'color','w');
if isfield(g,'lim')
    if ~isempty(g.lim)
        plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),line_style,'linewidth',2)
    end
end
contour(rPlot,zPlot,psiPlot,psiLevels,'-','linewidth',1);
contour(rPlot,zPlot,psiPlot,sepVal*[1,1],'k-','linewidth',2);
clim(sort([psi_min,psi_max]));
hc = colorbar; set(hc,'fontsize',14); colormap(turbo); hc.Label.String = cLabel;
xlabel('R (m)','fontsize',14); ylabel('Z (m)','fontsize',14)
axis equal; axis tight;
ax = axis;
if isfield(g,'gfilename')
    h=text(ax(1)+(ax(2)-ax(1))*0.48,ax(3)+(ax(4)-ax(3))*0.97,g.gfilename,'fontsize',8);
    set(h,'interpreter','none');
end

plot(g.rmaxis,g.zmaxis,'bo')
% Find the xpt(s)
haveXpt = 0;
if isfield(g,'bdry') && ~isempty(g.bdry) && any(g.bdry(1,:) > 1e-4)
    xpt_info = find_xpt_jl(g,1,1,1e-8,1);
    xr1 = xpt_info.rx;
    xz1 = xpt_info.zx;
    xr2 = xpt_info.rx2;
    xz2 = xpt_info.zx2;
    switch lower(plot_psi)
        case lower('psiN')
            psi_this = calc_psiN(g,xr2,xz2);           
        case lower('psi')
            psi_this = calc_psi(g,xr2,xz2);            
    end
    if ~isnan(psi_this)
        haveXpt = 1;
        contour(rPlot,zPlot,psiPlot,[1,1]*psi_this,'k-','linewidth',2)
        plot(xr1,xz1,'bx'); text(xr1+0.01,xz1,'x1','fontsize',8)
        plot(xr2,xz2,'b*'); text(xr2+0.02,xz2,'x2','fontsize',8)
    end
end

if plot_b_components
    [rMesh,zMesh] = meshgrid(rPlot,zPlot);
    Bout = bfield_geq_bicub(g,rMesh(:).',zMesh(:).',quiet);
    Br = reshape(Bout.br,size(rMesh));
    Bz = reshape(Bout.bz,size(rMesh));
    Btor = reshape(Bout.bphi,size(rMesh));
    Bpol = sqrt(Br.^2 + Bz.^2);

    bPlot(:,:,1) = Br;
    bPlot(:,:,2) = Bpol;
    bPlot(:,:,3) = Btor;
    bLabels = {'B_R (T)','B_{pol} (T)','B_\phi (T)'};
    nBPlots = 3;
    if debug_bpol_est
        if isfield(g,'Bpol_est') && ~isempty(g.Bpol_est)
            Bpol_est = interp2(g.r,g.z,g.Bpol_est.',rMesh,zMesh);
        else
            if isempty(g.cpasma)
                Bpol_est = nan(size(rMesh));
            else
                mu0 = 4*pi*1e-7;
                rminor = sqrt((rMesh-g.rmaxis).^2 + (zMesh-g.zmaxis).^2);
                Bpol_est = mu0*abs(g.cpasma)./(2*pi*rminor);
                Bpol_est(rminor == 0) = nan;
            end
        end
        bPlot(:,:,4) = Bpol_est;
        bLabels{4} = 'B_{pol,est} (T)';
        nBPlots = 4;
        Bpol_clim = [min(Bpol(:),[],'omitnan'),max(Bpol(:),[],'omitnan')];
    end

    figure; set(gcf,'color','w');
    if debug_bpol_est
        tiledlayout(2,2);
    else
        tiledlayout(1,3);
    end
    for i = 1:nBPlots
        nexttile; hold on; box on;
        btmp = bPlot(:,:,i);
        imagesc(rPlot,zPlot,btmp); set(gca,'ydir','normal');
        if any(~isnan(btmp(:)))
            if i == 1
                bmax = max(abs(btmp(:)),[],'omitnan');
                clim([-bmax,bmax]);
            elseif i == 4
                clim(Bpol_clim);
            else
                clim([min(btmp(:),[],'omitnan'),max(btmp(:),[],'omitnan')]);
            end
        end
        contour(rPlot,zPlot,psiPlot,psiLevels,'k-','linewidth',0.5);
        contour(rPlot,zPlot,psiPlot,sepVal*[1,1],'k-','linewidth',2);
        if isfield(g,'lim')
            if ~isempty(g.lim)
                plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),line_style,'linewidth',2)
            end
        end
        if haveXpt
            contour(rPlot,zPlot,psiPlot,[1,1]*psi_this,'k-','linewidth',2)
            plot(xr1,xz1,'bx'); text(xr1+0.01,xz1,'x1','fontsize',8)
            plot(xr2,xz2,'b*'); text(xr2+0.02,xz2,'x2','fontsize',8)
        end
        plot(g.rmaxis,g.zmaxis,'bo')
        xlabel('R (m)','fontsize',14); ylabel('Z (m)','fontsize',14)
        title(bLabels{i},'fontsize',14)
        axis equal; axis tight;
        hc = colorbar; set(hc,'fontsize',14); hc.Label.String = bLabels{i};
    end
    colormap(turbo);
end


end

%%
function [r2,z2,p2] = refine_psi(fac,g,plot_psi)
r2 = linspace(g.r(1),g.r(end),length(g.r)*fac);
z2 = linspace(g.z(1),g.z(end),length(g.z)*fac);
for i = 1:length(z2)
    switch lower(plot_psi)
        case lower('psiN')
            p2(i,:) = calc_psiN(g,r2,z2(i)*ones(size(r2)));
        case lower('psi')
            p2(i,:) = calc_psi(g,r2,z2(i)*ones(size(r2)));
        otherwise
            error('Bad value for plot_psi')
    end
end
end
