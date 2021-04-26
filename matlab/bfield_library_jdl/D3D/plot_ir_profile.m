clearvars;
elm_file = [];

% fname = 'C:\Work\DIII-D\APS 2016\irtv_data\irtv_156855_2000_4900.mat'; shot = 156855; SHIFT_CM = 3.5;
% fname = 'C:\Work\DIII-D\APS 2016\irtv_data\irtv_156861_2000_4900.mat'; shot = 156861; SHIFT_CM = 3.5;
% fname = 'C:\Work\DIII-D\APS 2016\irtv_data\irtv_156867_2000_4900.mat'; shot = 156867; SHIFT_CM = 0;
% fname = 'C:\Work\DIII-D\APS 2016\irtv_data\irtv_156871_2000_4900.mat'; shot = 156871; SHIFT_CM = 0;
fname = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\irtv_174306_3000_4000.mat'; shot = 174306; SHIFT_CM = 1.5;



% twin = [4300,4400];
% twin = [3300,3700];
% twin = [3200,4800];
tWinMS = [3100,4000];
elmWin = [0.8,0.95];
% elm_win = [0.9,0.95];



sparseifyLines = 20;    % Step for time trace indices, only in plotting
crangeElmPlot = 2;  % Colorbar control for ELM cycle plot: 1 = [0,1], else tight to elmWin set above

elmThreshNorm = 0.3;   % In normalized magnitude for ELM filter
dtThreshMS = 1;  % Discard "ELM" spikes repeated in this threshold (ms)
% Select ir radius to find ELMS (in dR from Outer SP)
dROutFindELMCM = -10;

WRITE_FILE = 1;
MAKE2DPLOTS = 0;

dROuterWantRangeCM = [-10,20];
dRInnerWantRangeCM = [-5,8];
nOuterBinIR = 12;  % Number of radial bins in dR range
nInnerBinIR = 12;  % Number of radial bins in dR range

% fname_stub = 'out';



%-------------------------------------------------------------------------%

%% Load file(s)
IR = get_IR_data(fname,tWinMS);
IR = filterDataElms(IR,dROutFindELMCM,elm_file,elmThreshNorm,dtThreshMS);
if MAKE2DPLOTS
    plotIR2D(IR,dROuterWantRangeCM,dRInnerWantRangeCM);
end
IR = elmFilter(IR,elmWin);
% IR = binIR(IR,dROuterWantRangeCM,dRInnerWantRangeCM,nBinIR);
IR.Binned.Outer = binIR(IR.ElmFiltered.time,IR.ElmFiltered.dROutCM,IR.ElmFiltered.dataMWm2,dROuterWantRangeCM,nOuterBinIR);
IR.Binned.Inner = binIR(IR.ElmFiltered.time,IR.ElmFiltered.dRInCM,IR.ElmFiltered.dataMWm2,dRInnerWantRangeCM,nInnerBinIR);

%% Plot data
plotVsElmCycle(IR,IR.ElmFiltered.dROutCM,IR.Binned.Outer,crangeElmPlot,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot);
plotVsElmCycle(IR,IR.ElmFiltered.dRInCM ,IR.Binned.Inner,crangeElmPlot,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot);

plotVsTime(IR,IR.ElmFiltered.dROutCM,IR.Binned.Outer,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot);
plotVsTime(IR,IR.ElmFiltered.dRInCM ,IR.Binned.Inner,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot);


drawnow;

if WRITE_FILE
    [pathstr] = fileparts(fname);    
    fname_out = sprintf('IR_profile_%d_twin_%d_%d_elm_%d_%d',shot,floor(tWinMS(1)),floor(tWinMS(2)),floor(elmWin(1)*100),floor(elmWin(2)*100));
    fullFileNameOut = fullfile(pathstr,[fname_out,'.mat']);
    fprintf('Writing file: %s\n',fullFileNameOut)
    writeFile(IR,fullFileNameOut,SHIFT_CM);
end

function writeFile(fullIR,fullFileNameOut,SHIFT_CM)

% Write file



%     %% dat

%     fullFileNameOut = fullfile(pathstr,[fname_out,'.dat']);
%     fid = fopen(fullFileNameOut,'w');
%     fprintf('Writing file: %s\n',fullFileNameOut)
%     fprintf(fid,'%d\n',length(dRBinCensCM));
%     for i = 1:length(dRBinCensCM)
%         fprintf(fid,'%e %e %e\n',dRBinCensCM(i),IR.Binned.Outer.dataMWm2(i),IR.Binned.Outer.errMWm2(i));
%     end
%     fclose(fid);

%% mat



IR.Outer.dRSep = fullIR.Binned.Outer.dRBinCensCM./100;
IR.Outer.dRSepShift = SHIFT_CM/100;
IR.Outer.q = fullIR.Binned.Outer.dataMWm2*1e6;
IR.Outer.qErr = fullIR.Binned.Outer.errMWm2*1e6;

IR.Inner.dRSep = fullIR.Binned.Inner.dRBinCensCM./100;
IR.Inner.dRSepShift = SHIFT_CM/100;
IR.Inner.q = fullIR.Binned.Inner.dataMWm2*1e6;
IR.Inner.qErr = fullIR.Binned.Inner.errMWm2*1e6;

save(fullFileNameOut,'IR')


end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function plotVsTime(IR,dRElmCM,IRBinned,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot)
% PLOT WITH COLOR = TIME

num_lines = size(IR.ElmFiltered.dataMWm2,2);
num_col = 1024;
mycol = colorflipper(num_col,'parula');

crange = tWinMS;
figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:sparseifyLines:num_lines
    icuse = interp1(linspace(crange(1),crange(2),num_col),1:num_col,IR.ElmFiltered.time(i));
    cuse = mycol(round(icuse),:);
    plot(dRElmCM(:,i)+SHIFT_CM,IR.ElmFiltered.dataMWm2(:,i),'-','color',cuse)
end
colorbar;
set(gca,'clim',crange);
% errorbar(drout_bin,IR.Binned.Outer.dataMWm2,IR.Binned.Outer.errMWm2,'ko--','linewidth',2)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2,'k','linewidth',3)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2-IRBinned.errMWm2,'k--','linewidth',2)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2+IRBinned.errMWm2,'k--','linewidth',2)
xlim = [IRBinned.dRBinEdgesCM(1)*0.99,IRBinned.dRBinEdgesCM(end)*1.01];
set(gca,'xlim',xlim)
ylim = get(gca,'ylim');
set(gca,'ylim',[0,ylim(2)]);
xlabel('R-R_{sep}^{out} [cm]')
ylabel('Q [MW/m^2]')
title(sprintf('Shot = %d, twin = [%3.1f,%3.1f], elm = [%d,%d] percent',shot,tWinMS(1)/1000,tWinMS(2)/1000,floor(elmWin(1)*100),floor(elmWin(2)*100)))
text(xlim(1)+diff(xlim)*0.6,ylim(2)*0.95,sprintf('Ncurve=%d',size(IR.ElmFiltered.dataMWm2,2)),'fontsize',14)

end

function plotVsElmCycle(IR,dRElmCM,IRBinned,crangeElmPlot,elmWin,tWinMS,sparseifyLines,SHIFT_CM,shot)

num_lines = size(IR.ElmFiltered.dataMWm2,2);
num_col = 1024;
mycol = colorflipper(num_col,'parula');

if crangeElmPlot == 1
    crange = [0,1];
else
    crange = elmWin;
end

figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:sparseifyLines:num_lines
    icuse = interp1(linspace(crange(1),crange(2),num_col),1:num_col,IR.ElmFiltered.elmCycle(i));
    cuse = mycol(round(icuse),:);
    plot(dRElmCM(:,i)+SHIFT_CM,IR.ElmFiltered.dataMWm2(:,i),'-','color',cuse)
end
colorbar;
set(gca,'clim',crange);
% errorbar(drout_bin,IR.Binned.Outer.dataMWm2,IR.Binned.Outer.errMWm2,'ko--','linewidth',2)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2,'k','linewidth',3)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2-IRBinned.errMWm2,'k--','linewidth',2)
plot(IRBinned.dRBinCensCM+SHIFT_CM,IRBinned.dataMWm2+IRBinned.errMWm2,'k--','linewidth',2)
xlim = [IRBinned.dRBinEdgesCM(1)*0.99,IRBinned.dRBinEdgesCM(end)*1.01];
set(gca,'xlim',xlim)
ylim = get(gca,'ylim');
set(gca,'ylim',[0,ylim(2)]);
xlabel('R-R_{sep}^{out} [cm]')
ylabel('Q [MW/m^2]')
title(sprintf('Shot = %d, twin = [%3.1f,%3.1f], elm = [%d,%d] percent',shot,tWinMS(1)/1000,tWinMS(2)/1000,floor(elmWin(1)*100),floor(elmWin(2)*100)))
text(xlim(1)+diff(xlim)*0.6,ylim(2)*0.95,sprintf('Ncurve=%d',size(IR.ElmFiltered.dataMWm2,2)),'fontsize',14)
end


function IR = get_IR_data(fname,tWinMS)
myir = load(fname);
ir = myir.ir;
% RAW
timeMS      = double(ir.time);
data        = double(ir.data);
drout       = double(ir.drout);
drin        = double(ir.drin);
IR.radiusCM = double(ir.radius);

% TIME SELECT
t_good = find(timeMS >= tWinMS(1) & timeMS <= tWinMS(2));
fprintf('Found %d time points in time window [%d,%d] (ms)\n',length(t_good),tWinMS(1),tWinMS(2));

IR.timeMS = timeMS(t_good);
IR.nTimes = length(timeMS(t_good));
IR.dataMWm2 = data(:,t_good);
IR.dROutCM = drout(:,t_good);
IR.dRInCM = drin(:,t_good);

end



function IR = filterDataElms(IR,dROutFindELMCM,elm_file,elmThreshNorm,dtThreshMS)

%% ELM Filter
if ~isempty(elm_file)
    % Provide something like FS data
    % Careful, haven't checked this recently
    error('confirm this works')
    elm = check_elm_timing(elm_file,tWinMS);
    IR.TimeELMCycle = -1*ones(size(IR.timeMS));
    for i=1:length(IR.timeMS)
        ind  = find(elm.telm(1:end-1) <= IR.timeMS(i) & elm.telm(2:end) > IR.timeMS(i));
        if ~isempty(ind)
            dte = elm.telm(ind+1)-elm.telm(ind);
            dts = IR.timeMS(i) - elm.telm(ind);
            IR.TimeELMCycle(i) = dts/dte;
        end
    end
else
    %% Use IR data at a given radius
    dat_1d = zeros(size(IR.timeMS));
    
    % Find the index of the ELM filter radius
    for i = 1:length(IR.timeMS)
        [~,ix] = min(abs(IR.dROutCM(:,i) - dROutFindELMCM));
        dat_1d(i) = IR.dataMWm2(ix,i);
    end
    % Idenfity spikes and spike times
    MedFiltIR = 1;
    Spikes = find_spikes(IR.timeMS,dat_1d,elmThreshNorm,dtThreshMS,MedFiltIR);
    % Assign portion of ELM cycle to each IR time
    IR.TimeELMCycle = -1*ones(size(IR.timeMS));
    for i=1:length(IR.timeMS)
        ind2 = find(Spikes.times(1:end-1) <= IR.timeMS(i) & Spikes.times(2:end) > IR.timeMS(i));
        if ~isempty(ind2)
            dte2 = Spikes.times(ind2+1)-Spikes.times(ind2);
            dts2 = IR.timeMS(i) - Spikes.times(ind2);
            IR.TimeELMCycle(i) = dts2/dte2;
        end
    end
end
end

function plotIR2D(IR,dROuterWantRangeCM,dRInnerWantRangeCM)
%
% IR PLOT 2D
%

figure; hold on; box on;
surf(IR.timeMS,IR.dROutCM,IR.dataMWm2,'edgecolor','none')
colorbar;
set(gca,'ylim',[dROuterWantRangeCM(1)*0.99,dROuterWantRangeCM(2)*1.01])
ylabel('R-R_{sep}^{out} [cm]')
title('Q vs dR outer [MW/m^2]')
xlabel('time [ms]')

figure; hold on; box on;
surf(IR.timeMS,IR.dRInCM,IR.dataMWm2,'edgecolor','none')
colorbar;
set(gca,'ylim',[dRInnerWantRangeCM(1)*0.99,dRInnerWantRangeCM(2)*1.01])
ylabel('R-R_{sep}^{in} [cm]')
title('Q vs dR inner [MW/m^2]')
xlabel('time [ms]')

% Plot vs R with Rsep included
figure; hold on; box on;
hp = pcolor(IR.timeMS,IR.radiusCM,IR.dataMWm2);
set(hp,'edgecolor','none');
RSepOut = IR.radiusCM(1)-IR.dROutCM(1,:);  % Indices don't matter here
RSepIn  = IR.radiusCM(1)+IR.dRInCM(1,:);
plot(IR.timeMS,RSepOut,'k')
plot(IR.timeMS,RSepIn,'k')
colorbar;
ylabel('R [cm]')
title('Q vs R [MW/m^2]')
xlabel('time [ms]')
end

function IR = elmFilter(IR,elmWin)

% ELM CYCLE SELECT
inElmCycle = find(IR.TimeELMCycle >= elmWin(1) & IR.TimeELMCycle <= elmWin(2));
IR.ElmFiltered.elmCycle = IR.TimeELMCycle(inElmCycle);
IR.ElmFiltered.time     = IR.timeMS(inElmCycle);
IR.ElmFiltered.dataMWm2 = IR.dataMWm2(:,inElmCycle);
IR.ElmFiltered.dROutCM  = IR.dROutCM(:,inElmCycle);
IR.ElmFiltered.dRInCM   = IR.dRInCM(:,inElmCycle);
end

function IRBinned = binIR(timeMS,radiusCM,dataMWm2,dRBinRange,nBin)

% Bins
IRBinned.dRBinEdgesCM = linspace(dRBinRange(1),dRBinRange(2),nBin+1);
IRBinned.dRBinCensCM  = (IRBinned.dRBinEdgesCM(2:end) +IRBinned.dRBinEdgesCM(1:end-1))/2;

IRBinned.dataMWm2 = nan(1,nBin);
IRBinned.errMWm2 = nan(1,nBin);
for ibin = 1:nBin
    databin = nan(1,length(timeMS));
    databinvar = nan(1,length(timeMS));
    for jtime = 1:length(timeMS)
        inBin = find( radiusCM(:,jtime) >= IRBinned.dRBinEdgesCM(ibin) & radiusCM(:,jtime) < IRBinned.dRBinEdgesCM(ibin+1));
        this = dataMWm2(inBin,jtime);
        databin(jtime)    = mean(this);
        databinvar(jtime) = var(this);
    end
    %     IRBinned.dataMWm2(ibin) = sum(databin./databinvar)./sum(1./databinvar);
    %     IRBinned.errMWm2(ibin) = sqrt(1./sum(1./databinvar));
    
    IRBinned.dataMWm2(ibin) = mean(databin);
    %     IR.Binned.Outer.errMWm2(ibin) = std(databin)./sqrt(length(databin));
    IRBinned.errMWm2(ibin) = std(databin);
    
    %     IRBinned.errMWm2(ibin) = std(databin);
end
end