clearvars;

% Channels 1:8 are on shelf, increasing in radius
% 9:19 are inner div, moving down center stack then out in radius

SAVE_IT = 0;
tWinMS = [3200,4000];
% elmWin = [0.0,1.0];
elmWin = [0.3,0.7];
crangeElmPlot = 2;  % Colorbar control for ELM cycle plot: 1 = [0,1], else tight to elmWin set above
elmThreshNorm = 0.1;
elmTimeNoRepeatMS = 1;
% signalNamePlot = 'heatflux';  % Just used to make figures
signalNamePlot = 'jsat';  % Just used to make figures
% signalNamePlot = 'dens';  % Just used to make figures
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\power13mw\g174310.03500_153';

fileNameLP = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\LP_174306.mat';
elmFile = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\fs_174306_2000_5000.mat';

% fileNameLP = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\LP_174310.mat';
% elmFile = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\fs_174310_2000_5000.mat';

[filepath,name,ext] = fileparts(fileNameLP);
outfile = fullfile(filepath,[name,'_processed',ext]);

%------------
chanShelf = 1:8;
chanInner = 9:19;

LP = load(fileNameLP); LP = LP.lp;
% g = readg_g3d(gfile_name);


%% Filter on tWin
LP = getLPDataTWin(LP,tWinMS);
LP = divideLPShelfInner(LP,chanShelf,chanInner);
LP = makeLPRadialProf(LP);

elmData = load(elmFile); tElm = double(elmData.fs.t00); sElm = double(elmData.fs.fs00);
iStartElmFilt = find(tElm > tWinMS(1),1,'first');
iEndElmFilt = find(tElm > tWinMS(2),1,'first');
Spikes = find_spikes(tElm(iStartElmFilt:iEndElmFilt),sElm(iStartElmFilt:iEndElmFilt),elmThreshNorm,elmTimeNoRepeatMS,1);

LP = filterLPElms(LP,Spikes,elmWin);

figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:length(LP.Shelf.chan)
    thisChan = LP.Shelf.chan{i};
%     plot(thisChan.delrsepout,thisChan.temp,'.')    
    plot(thisChan.dLSepOut,thisChan.(signalNamePlot),'x')    
    xlabel('dL_{sep} Outer')
    ylabel(signalNamePlot)
    title('Shelf')
end
plot(LP.Profiles.Shelf.dLSepOut,medflt1d_jl(LP.Profiles.Shelf.(signalNamePlot),100))
% plot(LP.ProfilesELMFiltered.Shelf.dLSepOut,LP.ProfilesELMFiltered.Shelf.(signalNamePlot),'ko')
plot(LP.ProfilesELMFiltered.Shelf.dLSepOut,medflt1d_jl(LP.ProfilesELMFiltered.Shelf.(signalNamePlot),100),'k','linew',2)






num_col = 1024;

if crangeElmPlot == 1
    crange = [0,1];
else
    crange = elmWin;
end
icuse = interp1(linspace(crange(1),crange(2),num_col),1:num_col,LP.Profiles.Shelf.periodELMCycle);

figure; hold on; box on; grid on; set(gcf,'color','w');
iBad = LP.Profiles.Shelf.periodELMCycle == -1;
scatter(LP.Profiles.Shelf.dLSepOut(~iBad),LP.Profiles.Shelf.(signalNamePlot)(~iBad),12,LP.Profiles.Shelf.periodELMCycle(~iBad),'filled')
colorbar
xlabel('dL_{sep} Outer')
ylabel(signalNamePlot)
title('Shelf')
colormap(colorflipper(num_col,'plasma'));
set(gca,'clim',crange)


figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:length(LP.Inner.chan)
    thisChan = LP.Inner.chan{i};
%     plot(thisChan.delrsepin,thisChan.temp,'.')    
    plot(thisChan.dLSepIn,thisChan.(signalNamePlot),'x')   
    xlabel('dL_{sep} Inner')
    ylabel(signalNamePlot)
    title('Inner')
end
plot(LP.ProfilesELMFiltered.Inner.dLSepIn,medflt1d_jl(LP.ProfilesELMFiltered.Inner.(signalNamePlot),100),'k','linew',2)

if SAVE_IT
save(outfile,'LP');
end

function LP = filterLPElms(LP,Spikes,elmWin)

locations = {'Shelf','Inner'};
for i = 1:length(locations)
    thisLocationName = locations{i};

    %% Tag time in ELM cycle
    thisTimeSeries = LP.Profiles.(thisLocationName).time;
    LP.Profiles.(thisLocationName).periodELMCycle = -1*ones(size(thisTimeSeries));
    
    for j=1:length(thisTimeSeries)
        ind2 = find(Spikes.times(1:end-1) <= thisTimeSeries(j) & Spikes.times(2:end) > thisTimeSeries(j));
        if ~isempty(ind2)
            dte2 = Spikes.times(ind2+1)-Spikes.times(ind2);
            dts2 = thisTimeSeries(j) - Spikes.times(ind2);
            LP.Profiles.(thisLocationName).periodELMCycle(j) = dts2/dte2;
        end
    end    
    
    % Filter using elmWin
    fieldNames = fieldnames(LP.Profiles.(thisLocationName));
    inElmCycle = find(LP.Profiles.(thisLocationName).periodELMCycle >= elmWin(1) & LP.Profiles.(thisLocationName).periodELMCycle <= elmWin(2));
    for k = 1:length(fieldNames)
        thisFieldName = fieldNames{k};        
        LP.ProfilesELMFiltered.(thisLocationName).(thisFieldName) = LP.Profiles.(thisLocationName).(thisFieldName)(inElmCycle);
    end
    
end

end

function LP = divideLPShelfInner(LP,chanShelf,chanInner)

iCount = 1;
for i = chanShelf
    LP.Shelf.chan{iCount} = LP.chan{i};
    LP.Shelf.chan{iCount}.probeNumber = i;    
    iCount = iCount + 1;
end

iCount = 1;
for i = chanInner
    LP.Inner.chan{iCount} = LP.chan{i};
    LP.Inner.chan{iCount}.probeNumber = i;      
    iCount = iCount + 1;
end
end


function LP = makeLPRadialProf(LP)
    
%For now use Shelf for outer SP and Inner for inner

%% Combine into one array
locations = {'Shelf','Inner'};
for j = 1:length(locations)
    thisLocation = locations{j};
    indexOffset = 0;
    for i = 1:length(LP.(thisLocation).chan)
        thisChan = LP.(thisLocation).chan{i};
        nD = thisChan.ntimes;
        theseInds = 1+indexOffset:indexOffset+nD;
        fieldNames = fieldnames(thisChan);
        for k = 1:length(fieldNames)
            thisFieldName = fieldNames{k};
            if length(thisChan.(thisFieldName)) == nD
                LP.Profiles.(thisLocation).(thisFieldName)(theseInds) = thisChan.(thisFieldName);
            end            
        end
        indexOffset = indexOffset + nD;
    end
        
end

%% sort by dLSep
[~,iSort] = sort(LP.Profiles.Shelf.dLSepOut);
fieldNames = fieldnames(LP.Profiles.Shelf);
for i = 1:length(fieldNames)
    LP.Profiles.Shelf.(fieldNames{i}) = LP.Profiles.Shelf.(fieldNames{i})(iSort);
end

[~,iSort] = sort(LP.Profiles.Inner.dLSepIn);
fieldNames = fieldnames(LP.Profiles.Inner);
for i = 1:length(fieldNames)
    LP.Profiles.Inner.(fieldNames{i}) = LP.Profiles.Inner.(fieldNames{i})(iSort);
end


end


function LPOut = getLPDataTWin(LP,tWinMS)

for i = 1:length(LP.chan)
    thisChan = LP.chan{i};
    if thisChan.ierr ~= 0
        error('bad channel?')
    end
    
    keepInds = find(thisChan.time >= tWinMS(1) & thisChan.time <= tWinMS(2));
    fieldNames = fieldnames(thisChan);
    for k = 1:length(fieldNames)
        thisFieldName = fieldNames{k};
        if length(thisChan.(thisFieldName)) == length(thisChan.time)
            LPOut.chan{i}.(thisFieldName) = double(thisChan.(thisFieldName)(keepInds));
        else
            LPOut.chan{i}.(thisFieldName) = thisChan.(thisFieldName);
        end        
    end
    LPOut.chan{i}.ntimes = length(keepInds);
    LPOut.chan{i}.dLSepOut = sqrt( LPOut.chan{i}.delrsepout.^2 + LPOut.chan{i}.delzsepout.^2);
    LPOut.chan{i}.dLSepIn  = sqrt( LPOut.chan{i}.delrsepin.^2  + LPOut.chan{i}.delzsepin.^2 );
end

end
