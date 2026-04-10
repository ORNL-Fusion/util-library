clearvars;


config_name = {'D3-6','D1-1','D2-2'};
config_name = {'D1-1'};

verbose = 1;
SIMPLIFY_COILS = 1;
PLOT_GEOMETRY = 0;
PLOT_COILS = 1;
PLOT_FLUX_CONTOURS = 0;
PLOT_LC = 'left';
PLOT_SCATTER = 1;
PLOT_TRIANGULATION = 1;
PLOT_LC_DEBUG = 0;
PLOT_LUFS = 0;
PROGRESS_OUTPUT = 1;
USE_PARFOR = 1;
TRIANGULATION_ALPHA = 0.08; % 0.05, 0.03, 0.08 --  bigger smoother - smaller tighter
nr_flux = 200;
nz_flux = 400;
Rmin_flux = 0;
Rmax_flux = 0.45;
ncontour_flux = 500;

dz_fieldline0 = 0.01;
USE_DL_FALLBACK = 1;
dl_fieldline0 = 0.01;
dl_length_factor = 1.5;
dcoil_cutoff_branch = 0.00001;
dl_max_point_jump = 0.2;


Geo = get_MPEX_geometry;

Zmin = min(Geo.Vessel.z);
Zmax = max(Geo.Vessel.z);

nlines_scattered = 350;
Rmin_scattered = 0.00;
Rmax_scattered = 0.15;
RSeed = linspace(Rmin_scattered,Rmax_scattered,nlines_scattered); 
ZSeed = Zmax*ones(size(RSeed));


Rmin_scattered = 0.00;
Rmax_scattered = 0.4;
RSeed2 = linspace(Rmin_scattered,Rmax_scattered,nlines_scattered); 
ZSeed2 = Zmin*ones(size(RSeed2));

RSeed = [RSeed,RSeed2];
ZSeed = [ZSeed,ZSeed2];

% ZSeed = [ZSeed,linspace(1,2.5,30)];
% RSeed = [RSeed,0.09*ones(1,30)];

% RSeed = [linspace(0,0.25,150), 0.10*ones(1,150)]; ZSeed = [Zmin*ones(1,150), linspace(Zmin,Zmax,150)];

%%
for i = 1:length(config_name)
    if PLOT_GEOMETRY
        figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)
    end

    %% Setup coils and bfield for this config
    if SIMPLIFY_COILS
        [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson_hybrid(config_name{i},[],3,3,verbose);
    else
        [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson(config_name{i},verbose);
    end
    Bfield.coil = Coil; Bfield.current = windingCurrent; Bfield.type = 'MPEX';



    %% Plot coils and vessel
    if PLOT_COILS || PLOT_GEOMETRY || PLOT_FLUX_CONTOURS
        [rcoil,zcoil] = get_coil_cross_sections(CoilGeometry);
    else
        rcoil = [];
        zcoil = [];
    end
    if PLOT_GEOMETRY
        if PLOT_COILS
            plot_coil_cross_section(rcoil,zcoil,0,currentPerWinding);
        end
        plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
        plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
        xlabel('Z (m)')
        ylabel('R (m)')
        title(sprintf('MPEX geometry: %s',config_name{i}))
        axis([-4,9,0,0.5])
    end

    if PLOT_FLUX_CONTOURS
        R1d_flux = linspace(Rmin_flux,Rmax_flux,nr_flux);
        Z1d_flux = linspace(Zmin,Zmax,nz_flux);
        [R2d_flux,Z2d_flux] = ndgrid(R1d_flux,Z1d_flux);
        psi2d = calc_psi_mpex(Coil,windingCurrent,R2d_flux,Z2d_flux);

        figure; set(gcf,'color','w')
        contour(Z2d_flux,R2d_flux,psi2d,ncontour_flux,'LineWidth',1)
        hold on
        plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
        plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
        if PLOT_COILS
            plot_coil_cross_section(rcoil,zcoil,0);
        end
        box on; grid on; set(gca,'fontsize',14)
        xlabel('Z (m)')
        ylabel('R (m)')
        title(sprintf('%s: Flux contours',config_name{i}))
        colorbar
        axis([-4,9,0,0.5])
    end

    %% Find LUFS
    if PLOT_LUFS
        f_lufs = find_lufs_MPEX(Bfield,Geo);
    else
        f_lufs = struct('r',[],'z',[],'hit_rz',[]);
    end


    %% Now let's compute connection length in space
    need_connection_plots = PLOT_SCATTER || PLOT_TRIANGULATION;


    need_left = strcmp(PLOT_LC,'left') || strcmp(PLOT_LC,'total');
    need_right = strcmp(PLOT_LC,'right') || strcmp(PLOT_LC,'total');
    need_total = strcmp(PLOT_LC,'total');
    
    if strcmp(PLOT_LC,'left')
        plot_title = sprintf('%s: Left connection length (m)',config_name{i});
    elseif strcmp(PLOT_LC,'right')
        plot_title = sprintf('%s: Right connection length (m)',config_name{i});
    elseif strcmp(PLOT_LC,'total')
        plot_title = sprintf('%s: Total connection length (m)',config_name{i});
    else
        error('Unknown PLOT_LC %s',PLOT_LC)
    end

    if need_connection_plots
        [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot,raw_fl_debug] = calc_connection_length_scattered(Bfield,Geo,RSeed,ZSeed,Zmin,Zmax,dz_fieldline0,need_left,need_right,need_total,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump,PROGRESS_OUTPUT,USE_PARFOR);
        
        if strcmp(PLOT_LC,'left')
            L_plot = L_backward_plot;
        elseif strcmp(PLOT_LC,'right')
            L_plot = L_forward_plot;
        else
            L_plot = L_total_plot;
        end
        
        if PLOT_SCATTER
            figure; set(gcf,'color','w')
            scatter(z_plot,r_plot,10,L_plot,'filled')
            hold on
            plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
            plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
            box on; grid on; set(gca,'fontsize',14)
            xlabel('Z (m)')
            ylabel('R (m)')
            title(plot_title)
            colorbar
            axis([-4,9,0,0.5])
            if PLOT_COILS
                plot_coil_cross_section(rcoil,zcoil,0);
            end
            if ~isempty(f_lufs.z)
                plot(f_lufs.z,f_lufs.r,'r')
            end
        end

        if PLOT_TRIANGULATION
            figure; set(gcf,'color','w')
            plot_triangulated_map(z_plot,r_plot,L_plot,Geo,plot_title,TRIANGULATION_ALPHA)
            if PLOT_COILS
                plot_coil_cross_section(rcoil,zcoil,0);
            end
            if ~isempty(f_lufs.z)
                plot(f_lufs.z,f_lufs.r,'r')
            end
        end

        if PLOT_LC_DEBUG
            figure; set(gcf,'color','w')
            scatter(z_plot,r_plot,10,L_plot,'filled')
            hold on
            for iseed = 1:numel(raw_fl_debug)
                if ~isempty(raw_fl_debug{iseed}.right.r)
                    plot(raw_fl_debug{iseed}.right.z,raw_fl_debug{iseed}.right.r,'k-','LineWidth',1)
                end
                if ~isempty(raw_fl_debug{iseed}.left.r)
                    plot(raw_fl_debug{iseed}.left.z,raw_fl_debug{iseed}.left.r,'k-','LineWidth',1)
                end
            end
            idl = method_plot == 2;
            if any(idl)
                scatter(z_plot(idl),r_plot(idl),28,'ko')
            end
            inonseed = seed_section_plot == 0;
            if any(inonseed)
                scatter(z_plot(inonseed),r_plot(inonseed),28,'mx')
            end
            plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
            plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
            box on; grid on; set(gca,'fontsize',14)
            xlabel('Z (m)')
            ylabel('R (m)')
            title(sprintf('%s [debug: circles=dl, x=non-seed vessel section]',plot_title))
            colorbar
            axis([-4,9,0,0.5])
            if PLOT_COILS
                plot_coil_cross_section(rcoil,zcoil,0);
            end
            if ~isempty(f_lufs.z)
                plot(f_lufs.z,f_lufs.r,'r')
            end
        end
    end

end

function f = find_lufs_MPEX(Bfield,Geo,RMax)

% if nargin < 3
    num_lines = 10;
% end
% if nargin < 5
    nbisect = 8;
% end

ZMin = min(Geo.Vessel.z) + 0.01;  % Below this Z lines are not checked for cutoff
if nargin < 3
    RMax = max(Geo.Target.r) + 0.01;
end
RMin = 1e-3;

% set up fieldline following
L = Geo.Target.z(1) - ZMin;
dz_fieldline = -0.01; 
nsteps = round(abs(L/dz_fieldline)); 
dz_fieldline = sign(dz_fieldline)*L/nsteps;

rr = linspace(RMin,RMax,num_lines);
phistart = zeros(size(rr));
fl = follow_fieldlines_rzphi_dz(Bfield,rr,Geo.Target.z(1),phistart,dz_fieldline,nsteps);
fl = clip_fl_at_vessel(fl,Geo.Vessel.r,Geo.Vessel.z);
isgood = ~isnan(fl.z(end,:));
ind_last_good_line = find(isgood,1,'last');

if isempty(ind_last_good_line)
    error('Could not identify LUFS bracket: all coarse lines hit the vessel before reaching ZMin.')
end
if ind_last_good_line == numel(rr)
    error('Could not identify LUFS bracket: all coarse lines passed through. Increase RMax.')
end

r_good = rr(ind_last_good_line);
r_bad = rr(ind_last_good_line+1);
fl_good.r = fl.r(:,ind_last_good_line);
fl_good.z = fl.z(:,ind_last_good_line);
fl_good.phi = fl.phi(:,ind_last_good_line);
fl_bad.r = fl.r(:,ind_last_good_line+1);
fl_bad.z = fl.z(:,ind_last_good_line+1);
fl_bad.phi = fl.phi(:,ind_last_good_line+1);

for ib = 1:nbisect
    r_mid = 0.5*(r_good + r_bad);
    fl_mid = trace_lufs_line(Bfield,Geo,r_mid,Geo.Target.z(1),dz_fieldline,nsteps);
    if ~isnan(fl_mid.z(end))
        r_good = r_mid;
        fl_good = fl_mid;
    else
        r_bad = r_mid;
        fl_bad = fl_mid;
    end
end

ind_hit = find(~isnan(fl_bad.z),1,'last');
if isempty(ind_hit)
    hit_rz = [NaN,NaN];
else
    hit_rz = [fl_bad.r(ind_hit),fl_bad.z(ind_hit)];
end

f.r = fl_good.r;
f.z = fl_good.z;
f.hit_rz = hit_rz;


end


function fl_single = trace_lufs_line(Bfield,Geo,r0,z0,dz_fieldline,nsteps)
fl_raw = follow_fieldlines_rzphi_dz(Bfield,r0,z0,0,dz_fieldline,nsteps);
fl_single = clip_fl_at_vessel(fl_raw,Geo.Vessel.r,Geo.Vessel.z);
end



function fout = clip_fl_at_vessel(f,vessel_r,vessel_z)
% Clip field lines at the first point that leaves the vessel
% Same field line structure returned, but points after vessel exit are set
% to NaN
fout.r = NaN(size(f.r));
fout.z = NaN(size(f.z));
fout.phi = NaN(size(f.phi));

for i = 1:size(f.r,2)
    isin2 = inpolygon(f.r(:,i),f.z(:,i),vessel_r,vessel_z);
    % Keep points up to the last point still inside the vessel
    is2 = find(isin2 == 0,1,'first') - 1;
    if isempty(is2)
        is2 = length(isin2);
    end
    if is2 > 0
        fout.r(1:is2,i) = f.r(1:is2,i);
        fout.z(1:is2,i) = f.z(1:is2,i);
        fout.phi(1:is2,i) = f.phi(1:is2,i);
    end
end
end


function fout = mask_fl_at_vessel(f,vessel_r,vessel_z)
% Same field line structure returned, but points outside the vessel are NaN
fout = f;
isin = inpolygon(f.r,f.z,vessel_r,vessel_z);
fout.r(~isin) = NaN;
fout.z(~isin) = NaN;
if isfield(f,'phi')
    fout.phi(~isin) = NaN;
end
end


function [Lc_forward,Lc_backward,Lc_total] = calc_connection_length_sections(f,calc_forward,calc_backward,calc_total)
% Compute connection length within each continuous in-vessel section.
if nargin < 2
    calc_forward = 1;
end
if nargin < 3
    calc_backward = 1;
end
if nargin < 4
    calc_total = 1;
end

Lc_forward = NaN(size(f.r));
Lc_backward = NaN(size(f.r));
Lc_total = NaN(size(f.r));

for i = 1:size(f.r,2)
    isgood = ~isnan(f.r(:,i)) & ~isnan(f.z(:,i));
    if ~any(isgood)
        continue
    end

    idiff = diff([0;isgood;0]);
    istart = find(idiff == 1);
    iend = find(idiff == -1) - 1;
    for j = 1:length(istart)
        ind = istart(j):iend(j);
        Lseg = [0;cumsum(sqrt(diff(f.r(ind,i)).^2 + diff(f.z(ind,i)).^2))];
        if calc_backward
            Lc_backward(ind,i) = Lseg;
        end
        if calc_forward
            Lc_forward(ind,i) = Lseg(end) - Lseg;
        end
        if calc_total
            Lc_total(ind,i) = Lseg(end);
        end
    end
end
end


function [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot,raw_fl_debug] = calc_connection_length_scattered(Bfield,Geo,RSeed,ZSeed,Zmin,Zmax,dz_fieldline0,need_left,need_right,need_total,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump,progress_output,use_parfor)
if ~isequal(size(RSeed),size(ZSeed))
    error('RSeed and ZSeed must have the same size.')
end

if nargin < 16 || isempty(progress_output)
    progress_output = 0;
end
if nargin < 17 || isempty(use_parfor)
    use_parfor = 0;
end

nseed = numel(RSeed);
r_plot_cell = cell(nseed,1);
z_plot_cell = cell(nseed,1);
L_forward_plot_cell = cell(nseed,1);
L_backward_plot_cell = cell(nseed,1);
L_total_plot_cell = cell(nseed,1);
method_plot_cell = cell(nseed,1);
seed_section_plot_cell = cell(nseed,1);
raw_fl_debug = cell(nseed,1);

ztol = max(1e-9,1e-6*max(1,abs(Zmax - Zmin)));

if use_parfor
    if progress_output
        fprintf('Running %d connection-length lines with parfor\n',nseed);
        progress_queue = parallel.pool.DataQueue;
        progress_count = 0;
        afterEach(progress_queue,@update_parfor_progress);
    else
        progress_queue = [];
    end
    parfor iseed = 1:nseed
        [r_plot_cell{iseed},z_plot_cell{iseed},L_forward_plot_cell{iseed},L_backward_plot_cell{iseed}, ...
            L_total_plot_cell{iseed},method_plot_cell{iseed},seed_section_plot_cell{iseed},raw_fl_debug{iseed}] = ...
            calc_connection_length_single_seed(Bfield,Geo,RSeed(iseed),ZSeed(iseed),Zmin,Zmax,ztol,dz_fieldline0, ...
            need_left,need_right,need_total,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump);
        if progress_output
            send(progress_queue,iseed);
        end
    end
else
    for iseed = 1:nseed
        if progress_output
            fprintf('Working on line %d of %d\n',iseed,nseed);
        end
        [r_plot_cell{iseed},z_plot_cell{iseed},L_forward_plot_cell{iseed},L_backward_plot_cell{iseed}, ...
            L_total_plot_cell{iseed},method_plot_cell{iseed},seed_section_plot_cell{iseed},raw_fl_debug{iseed}] = ...
            calc_connection_length_single_seed(Bfield,Geo,RSeed(iseed),ZSeed(iseed),Zmin,Zmax,ztol,dz_fieldline0, ...
            need_left,need_right,need_total,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump);
    end
end

r_plot = vertcat(r_plot_cell{:});
z_plot = vertcat(z_plot_cell{:});
L_forward_plot = vertcat(L_forward_plot_cell{:});
L_backward_plot = vertcat(L_backward_plot_cell{:});
L_total_plot = vertcat(L_total_plot_cell{:});
method_plot = vertcat(method_plot_cell{:});
seed_section_plot = vertcat(seed_section_plot_cell{:});

    function update_parfor_progress(~)
        progress_count = progress_count + 1;
        fprintf('Completed line %d of %d\n',progress_count,nseed);
    end
end


function [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot,raw_fl_debug] = ...
    calc_connection_length_single_seed(Bfield,Geo,r0,z0,Zmin,Zmax,ztol,dz_fieldline0,need_left,need_right,need_total,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump)
start_on_left = abs(z0 - Zmin) <= ztol;
start_on_right = abs(z0 - Zmax) <= ztol;

if start_on_left
    follow_right = need_left || need_right || need_total;
    follow_left = 0;
elseif start_on_right
    follow_left = need_left || need_right || need_total;
    follow_right = 0;
else
    follow_right = need_right || need_total;
    follow_left = need_left || need_total;
end

r_plot = [];
z_plot = [];
L_forward_plot = [];
L_backward_plot = [];
L_total_plot = [];
method_plot = [];
seed_section_plot = [];

fl_right = [];
fl_left = [];
fl_right_raw = empty_fl_struct;
fl_left_raw = empty_fl_struct;
right_complete = 0;
left_complete = 0;
if follow_right
    [fl_right,fl_right_raw,right_complete,right_method] = trace_connection_branch(Bfield,Geo,r0,z0,Zmax,dz_fieldline0,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump);
    if right_complete
        [Lc_forward_right,Lc_backward_right,~] = calc_connection_length_sections(fl_right,1,1,0);
    else
        Lc_forward_right = [];
        Lc_backward_right = [];
    end
else
    Lc_forward_right = [];
    Lc_backward_right = [];
end

if follow_left
    [fl_left,fl_left_raw,left_complete,left_method] = trace_connection_branch(Bfield,Geo,r0,z0,Zmin,dz_fieldline0,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump);
    if left_complete
        [Lc_forward_left,Lc_backward_left,~] = calc_connection_length_sections(fl_left,1,1,0);
    else
        Lc_forward_left = [];
        Lc_backward_left = [];
    end
else
    Lc_forward_left = [];
    Lc_backward_left = [];
end

raw_fl_debug = struct('right',fl_right_raw,'left',fl_left_raw);

if follow_right && right_complete
    [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot] = append_branch_samples( ...
        r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot, ...
        fl_right,Lc_forward_right,Lc_backward_right,'right',need_left,need_right,need_total,0,right_method);
end

if follow_left && left_complete
    [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot] = append_branch_samples( ...
        r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot, ...
        fl_left,Lc_forward_left,Lc_backward_left,'left',need_left,need_right,need_total,follow_right,left_method);
end
end


function [fl_branch,fl_raw,branch_complete,method_used] = trace_connection_branch(Bfield,Geo,r0,z0,zend,dz_fieldline0,USE_DL_FALLBACK,dl_fieldline0,dl_length_factor,dcoil_cutoff_branch,dl_max_point_jump)
L = zend - z0;
nsteps = max(1,ceil(abs(L/dz_fieldline0)));
dz_fieldline = L/nsteps;
bfield_branch = Bfield;
bfield_branch.dcoil_cutoff = dcoil_cutoff_branch;
bfield_branch.nsym = 1;
branch_complete = 0;
method_used = 1;
try
    [fl_raw,ierr_dz] = follow_fieldlines_rzphi_dz(bfield_branch,r0,z0,0,dz_fieldline,nsteps,1);
    branch_complete = (ierr_dz == 0) && ~has_large_fl_jump(fl_raw,dl_max_point_jump);
catch ME
    if ~USE_DL_FALLBACK
        rethrow(ME)
    else
        ierr_dz = 1;
    end
end

if ~branch_complete && USE_DL_FALLBACK
    lmax = dl_length_factor*abs(zend - z0);
    nsteps_dl = max(1,ceil(lmax/dl_fieldline0));
    dl = sign(zend - z0)*lmax/nsteps_dl;
    [fl_try,ierr_dl] = follow_fieldlines_rzphi_dl(bfield_branch,r0,z0,0,dl,nsteps_dl,1);
    method_used = 2;
    if ierr_dl == 0
        if sign(zend - z0) > 0
            ireach = find(fl_try.z >= zend,1,'first');
        else
            ireach = find(fl_try.z <= zend,1,'first');
        end

        if ~isempty(ireach)
            if ireach == 1
                fl_raw.r = fl_try.r(1,:);
                fl_raw.z = fl_try.z(1,:);
                fl_raw.phi = fl_try.phi(1,:);
            else
                z1 = fl_try.z(ireach-1);
                z2 = fl_try.z(ireach);
                frac = (zend - z1)/(z2 - z1);
                frac = max(0,min(1,frac));
                r_end = fl_try.r(ireach-1,:) + frac*(fl_try.r(ireach,:) - fl_try.r(ireach-1,:));
                phi_end = fl_try.phi(ireach-1,:) + frac*(fl_try.phi(ireach,:) - fl_try.phi(ireach-1,:));
                fl_raw.r = [fl_try.r(1:ireach-1,:);r_end];
                fl_raw.z = [fl_try.z(1:ireach-1,:);zend];
                fl_raw.phi = [fl_try.phi(1:ireach-1,:);phi_end];
            end
            branch_complete = ~has_large_fl_jump(fl_raw,dl_max_point_jump);
        else
            fl_raw = fl_try;
            branch_complete = 0;
        end
    else
        fl_raw = fl_try;
        branch_complete = 0;
    end
end

fl_branch = mask_fl_at_vessel(fl_raw,Geo.Vessel.r,Geo.Vessel.z);
end


function value = first_finite_scalar(arr)
ivalid = find(isfinite(arr),1,'first');
if isempty(ivalid)
    value = NaN;
else
    value = arr(ivalid);
end
end


function fl = empty_fl_struct
fl = struct('r',[],'z',[],'phi',[]);
end


function tf = has_large_fl_jump(fl,jump_threshold)
if isempty(fl.r) || numel(fl.r) < 2 || isempty(jump_threshold) || jump_threshold <= 0
    tf = false;
    return
end

dr = diff(fl.r(:));
dz = diff(fl.z(:));
step_jump = sqrt(dr.^2 + dz.^2);
tf = any(step_jump > jump_threshold);
end


function [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot] = append_branch_samples( ...
    r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot,method_plot,seed_section_plot, ...
    fl_branch,Lbranch_forward,Lbranch_backward,branch_dir,need_left,need_right,need_total,omit_seed,method_used)
ivalid = find(~isnan(fl_branch.r) & ~isnan(fl_branch.z));
if isempty(ivalid)
    return
end
if omit_seed && numel(ivalid) > 1
    ivalid = ivalid(2:end);
end

rline = fl_branch.r(ivalid);
zline = fl_branch.z(ivalid);
Lforward_branch = Lbranch_forward(ivalid);
Lbackward_branch = Lbranch_backward(ivalid);
seed_connected_valid = false(size(ivalid));
isgood = ~isnan(fl_branch.r) & ~isnan(fl_branch.z);
idiff = diff([0;isgood;0]);
istart = find(idiff == 1);
iend = find(idiff == -1) - 1;
if ~isempty(istart)
    seed_connected_valid = ivalid >= istart(1) & ivalid <= iend(1);
end

if strcmp(branch_dir,'right')
    if need_right
        Lright = Lforward_branch;
    else
        Lright = NaN(size(Lforward_branch));
    end
    if need_left || need_total
        Lleft = Lbackward_branch;
    else
        Lleft = NaN(size(Lforward_branch));
    end
else
    if need_left
        Lleft = Lforward_branch;
    else
        Lleft = NaN(size(Lforward_branch));
    end
    if need_right || need_total
        Lright = Lbackward_branch;
    else
        Lright = NaN(size(Lforward_branch));
    end
end
if need_total
    Ltotal = Lleft + Lright;
else
    Ltotal = NaN(size(Lforward_branch));
end

if need_right
    Lforward_use = Lright;
else
    Lforward_use = NaN(size(Lright));
end
if need_left
    Lbackward_use = Lleft;
else
    Lbackward_use = NaN(size(Lleft));
end
if need_total
    Ltotal_use = Ltotal;
else
    Ltotal_use = NaN(size(Ltotal));
end

ikeep = isfinite(Lforward_use) | isfinite(Lbackward_use) | isfinite(Ltotal_use);
if ~any(ikeep)
    return
end

r_plot = [r_plot;rline(ikeep)];
z_plot = [z_plot;zline(ikeep)];
L_forward_plot = [L_forward_plot;Lforward_use(ikeep)];
L_backward_plot = [L_backward_plot;Lbackward_use(ikeep)];
L_total_plot = [L_total_plot;Ltotal_use(ikeep)];
method_plot = [method_plot;method_used*ones(nnz(ikeep),1)];
seed_section_plot = [seed_section_plot;seed_connected_valid(ikeep)];
end


function plot_triangulated_map(z_plot,r_plot,L_plot,Geo,plot_title_str,triangulation_alpha)
z_plot = z_plot(:);
r_plot = r_plot(:);
L_plot = L_plot(:);

dt = delaunayTriangulation(z_plot,r_plot);
tri = dt.ConnectivityList;
hold on
ztri = reshape(dt.Points(tri,1),size(tri));
rtri = reshape(dt.Points(tri,2),size(tri));
ztri_centroid = mean(ztri,2);
rtri_centroid = mean(rtri,2);

if ~isempty(triangulation_alpha)
    shp = alphaShape(z_plot,r_plot,triangulation_alpha);
    inside_alpha = inShape(shp,ztri_centroid,rtri_centroid);
else
    shp = alphaShape(z_plot,r_plot);
    shp.Alpha = criticalAlpha(shp,'one-region');
    inside_alpha = inShape(shp,ztri_centroid,rtri_centroid);
end

inside_vessel = inpolygon(ztri(:,1),rtri(:,1),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri(:,2),rtri(:,2),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri(:,3),rtri(:,3),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri_centroid,rtri_centroid,Geo.Vessel.z,Geo.Vessel.r);
inside_tri = inside_vessel & inside_alpha;
tri = tri(inside_tri,:);

patch('Faces',tri, ...
    'Vertices',dt.Points, ...
    'FaceVertexCData',L_plot, ...
    'FaceColor','interp', ...
    'EdgeColor','none');
plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
box on; grid on; set(gca,'fontsize',14)
xlabel('Z (m)')
ylabel('R (m)')
title(plot_title_str)
colorbar
axis([-4,9,0,0.5])
end

