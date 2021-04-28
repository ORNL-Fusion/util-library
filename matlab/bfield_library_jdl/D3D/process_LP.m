clearvars;

% Channels 1:8 are on shelf, increasing in radius
% 9:19 are inner div, moving down center stack then out in radius

tWinMS = [3100,4000];
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\power13mw\g174310.03500_153';
fileNameLP = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\LP_174306.mat';

[filepath,name,ext] = fileparts(fileNameLP);
outfile = fullfile(filepath,[name,'_processed',ext]);

%------------
chanShelf = 1:8;
chanInner = 9:19;

LP = load(fileNameLP); LP = LP.lp;
g = readg_g3d(gfile_name);


%% Filter on tWin
LP = getLPDataTWin(LP,tWinMS);
LP = divideLPShelfInner(LP,chanShelf,chanInner);
LP = makeLPRadialProf(LP);


figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:length(LP.Shelf.chan)
    thisChan = LP.Shelf.chan{i};
%     plot(thisChan.delrsepout,thisChan.temp,'.')    
    plot(thisChan.dLSepOut,thisChan.jsat,'x')    
    xlabel('dL_{sep} Outer')
    ylabel('j_{sat} ()')
    title('Shelf')
end

figure; hold on; box on; grid on; set(gcf,'color','w');
for i = 1:length(LP.Inner.chan)
    thisChan = LP.Inner.chan{i};
%     plot(thisChan.delrsepin,thisChan.temp,'.')    
    plot(thisChan.dLSepIn,thisChan.jsat,'x')   
    xlabel('dL_{sep} Inner')
    ylabel('j_{sat} ()')
    title('Inner')
end

save(outfile,'LP');


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
            indexOffset = indexOffset + nD;
        end
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
        LPOut.chan{i} = thisChan;
        LPOut.chan{i}.ntimes = length(keepInds);
        LPOut.chan{i}.time = double(thisChan.time(keepInds));
        LPOut.chan{i}.temp = double(thisChan.temp(keepInds));
        LPOut.chan{i}.temp_err = double(thisChan.temp_err(keepInds));
        LPOut.chan{i}.dens = double(thisChan.dens(keepInds));
        LPOut.chan{i}.dens_err = double(thisChan.dens_err(keepInds));
        LPOut.chan{i}.pot = double(thisChan.pot(keepInds));
        LPOut.chan{i}.pot_err = double(thisChan.pot_err(keepInds));
        LPOut.chan{i}.isat = double(thisChan.isat(keepInds));
        LPOut.chan{i}.jsat = double(thisChan.jsat(keepInds));
        LPOut.chan{i}.angle = double(thisChan.angle(keepInds));
        LPOut.chan{i}.area = double(thisChan.area(keepInds));
        LPOut.chan{i}.psin = double(thisChan.psin(keepInds));
        LPOut.chan{i}.delrsepin = double(thisChan.delrsepin(keepInds));
        LPOut.chan{i}.delrsepout = double(thisChan.delrsepout(keepInds));
        LPOut.chan{i}.delzsepin = double(thisChan.delzsepin(keepInds));
        LPOut.chan{i}.delzsepout = double(thisChan.delzsepout(keepInds));
        LPOut.chan{i}.csq = double(thisChan.csq(keepInds));
        LPOut.chan{i}.res_err = double(thisChan.res_err(keepInds));
        LPOut.chan{i}.heatflux = double(thisChan.heatflux(keepInds));
        
        LPOut.chan{i}.dLSepOut = sqrt( LPOut.chan{i}.delrsepout.^2 + LPOut.chan{i}.delzsepout.^2);
        LPOut.chan{i}.dLSepIn  = sqrt( LPOut.chan{i}.delrsepin.^2 + LPOut.chan{i}.delzsepin.^2);
    end 

end

% plot(LP.hf_r_rsep+Rsep,medflt1d_jl(LP.hf,5),'rx','HandleVisibility','off')