function irtv = plot_d3d_irtv_data(fname,elm_cycle_min,elm_cycle_max,plotit,quiet)
if nargin == 0
    fname = 'C:\Users\jjl\Dropbox (ORNL)\Elijah\165908\165908\irtv_165908_3000_3200.mat';
    elm_cycle_min = 0.5;
    elm_cycle_max = 0.9;
    plotit = 1;
end

if nargin < 4
    plotit = 1;
end
if nargin < 5
    quiet = 0;
end


irtemp = load(fname);
ir = irtemp.ir;

nt = length(ir.time);
tmin = ir.time(1);
tmax = ir.time(end);

if ~quiet
    fprintf('time: %.2f < t < %.2f, nt = %d\n',tmin,tmax,nt)
end

% figure; hold on; box on;
% plot(ir.time,ir.es)
% title('elm cycle')

% ir.es = elm cycle
% ncount_temp = length(find(diff(ir.es) < -0.8));  % just make sure there are some elms
in_elm_cycle = ir.es >= elm_cycle_min & ir.es <= elm_cycle_max;
n_in_elm_cycle = sum(in_elm_cycle);
if n_in_elm_cycle <= 0
    error('no elms?')
end
if ~quiet
    fprintf('keeping %d of %d points after ELM cycle filter\n',n_in_elm_cycle,nt)
end


irtv.qfilt_static = ir.data(:,in_elm_cycle);
irtv.qmean_static = mean(irtv.qfilt_static,2);
irtv.qstd_static  = std(irtv.qfilt_static,[],2);
irtv.R     = ir.radius;
ir.elm_filt = [elm_cycle_min,elm_cycle_max];



q = ir.data(:,in_elm_cycle);
nsmooth = 1000;

dr_min_cut_out = -15; % cm -- trim data due to poor fit to dr_sep far from sep
dr_max_cut_out = 50; % cm -- trim data due to poor fit to dr_sep far from sep
dr_out = ir.drout(:,in_elm_cycle);
dr_smooth_out = linspace(min(dr_out(:)),max(dr_out(:)),nsmooth);
ir_smooth_all_out = nan(nsmooth,n_in_elm_cycle);
for i = 1:n_in_elm_cycle    
    x = dr_out(:,i);
    y = q(:,i);
    igood = find(x >= dr_min_cut_out & x <= dr_max_cut_out);
    ir_smooth_all_out(1:nsmooth,i) = interp1(x(igood),y(igood),dr_smooth_out);
end
ir_smooth_out = mean(ir_smooth_all_out,2,'omitnan');


irtv.ir_out = ir_smooth_out;
irtv.q_out = ir_smooth_out;

% dr_min_cut_in = -15; % cm -- trim data due to poor fit to dr_sep far from sep
% dr_max_cut_in = 50; % cm -- trim data due to poor fit to dr_sep far from sep
% dr_in = ir.drin(:,in_elm_cycle);
% dr_smooth_in = linspace(min(dr_in(:)),max(dr_in(:)),nsmooth);
% ir_smooth_all_in = nan(nsmooth,n_in_elm_cycle);
% for i = 1:n_in_elm_cycle    
%     x = dr_in(:,i);
%     y = q(:,i);
%     igood = find(x >= dr_min_cut_in & x <= dr_max_cut_in);
%     ir_smooth_all_in(1:nsmooth,i) = interp1(x(igood),y(igood),dr_smooth_in);
% end
% ir_smooth_in = mean(ir_smooth_all_in,2,'omitnan');

warning('in this case the inner mapped data is poor so I skipped!')



if plotit
    cols = lines;
    figure; hold on; box on; grid on;
    plot(ir.radius,ir.data,'color','k','linew',0.1)
    % figure; hold on; box on; grid on;
    plot(ir.radius,irtv.qfilt_static,'color',cols(1,:),'linew',2)
    plot(ir.radius,irtv.qmean_static,'color',cols(2,:),'linew',3)
    plot(ir.radius,irtv.qmean_static-irtv.qstd_static,'-.','color',cols(2,:),'linew',3)
    plot(ir.radius,irtv.qmean_static+irtv.qstd_static,'-.','color',cols(2,:),'linew',3)
    xlabel('R (m)')
    ylabel('q (MW/m^2)')
    title(sprintf('ELM average %.0f-%.0f%%',elm_cycle_min*100,elm_cycle_max*100))
end