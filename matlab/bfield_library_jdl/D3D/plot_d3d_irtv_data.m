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
nfound = sum(in_elm_cycle);
if nfound <= 0
    error('no elms?')
end
if ~quiet
    fprintf('keeping %d of %d points after ELM cycle filter\n',nfound,nt)
end
irtv.qall  = ir.data;
irtv.elmc  = ir.es;
irtv.qfilt = ir.data(:,in_elm_cycle);
irtv.drout = ir.drout(:,in_elm_cycle)./100;
irtv.drin  = ir.drin(:,in_elm_cycle)./100;
irtv.qmean = mean(irtv.qfilt,2);
irtv.qstd  = std(irtv.qfilt,[],2);
irtv.R     = ir.radius;
ir.elm_filt = [elm_cycle_min,elm_cycle_max];

if plotit
    cols = lines;
    figure; hold on; box on; grid on;
    plot(ir.radius,ir.data,'color','k','linew',0.1)
    % figure; hold on; box on; grid on;
    plot(ir.radius,irtv.qfilt,'color',cols(1,:),'linew',2)
    plot(ir.radius,irtv.qmean,'color',cols(2,:),'linew',3)
    plot(ir.radius,irtv.qmean-irtv.qstd,'-.','color',cols(2,:),'linew',3)
    plot(ir.radius,irtv.qmean+irtv.qstd,'-.','color',cols(2,:),'linew',3)
    xlabel('R (m)')
    ylabel('q (MW/m^2)')
    title(sprintf('ELM average %.0f-%.0f%%',elm_cycle_min*100,elm_cycle_max*100))
end