function [geo] = get_Proto_geometry(plotit,newfig,add_skimmer,target_position,add_sleeve)
% target_position = 1, target is at 7.5  (+1")
% target_position = 2, target is at 11.5 (+5")
if nargin < 1
    plotit = 0;
end
if nargin < 2
    newfig = 1;
end
% if nargin < 3 
%     add_skimmer = 0;
% end
% if nargin < 4
%     target_position = 1;
% end
% if nargin < 5
%     add_sleeve = 0;
% end

in2m = 0.0254;

% Get coil info for positioning elements
[rcoil,zcoil] = get_coil_cross_sections;
cmax = max(zcoil,[],2);  % max and min Z of the coils
cmin = min(zcoil,[],2);

% Possible target positions
if target_position == 1
    target.z = mean([cmax(7),cmin(8)]) + 1*in2m;  % midway between 7-8, 1" towards 8
    target.r = 3.75*in2m/2;  % diameter is 3.75 inches
elseif target_position == 2
    target.z = mean([cmax(11),cmin(12)]) + 5*in2m;  % 11.5 + 5"
    target.r = 4.5*in2m/2;  %4.5" diameter
else
    error(['Did not recognize target_position:',num2str(target_position)])
end

% "Sleeve"
sleeve.r = 0.08/2;  % ID = 80mm
sleeve.z1 = cmin(7);
sleeve.z2 = sleeve.z1 + 24*in2m; %upstream coil 7, to down 9.  24" 
sleeve.z = [sleeve.z1,sleeve.z2];

% Helicon
helicon_L  = 11.8*in2m;                    % L is 11.8" from Meitner
helicon_Z_cen = mean([cmax(3),cmin(4)]);   % Centered between coils 3 and 4
helicon.z1 = helicon_Z_cen - helicon_L/2;         
helicon.z2 = helicon_Z_cen + helicon_L/2;
helicon.zmid = helicon_Z_cen;
helicon.r = 4.95*in2m/2;                    % diameter is 5.1 inches --> update 4/21/16.  Nominal is ~4.95
helicon.z = [helicon.z1,helicon.z2];

% Skimmer
skimmer.z1 = (cmin(5) + cmax(5))/2 + 0.0702 - 0.5e-2;  % Finite width for finding intersections
skimmer.z2 = (cmin(5) + cmax(5))/2 + 0.0702 + 0.5e-2;  
skimmer.zmid = (skimmer.z1+skimmer.z2)/2;
skimmer.r = 0.058/2;
skimmer.z = [skimmer.z1,skimmer.z2];

dump.z = mean([cmax(11),cmin(12)]) - 3.70;  % 370cm from 11.5
dump.OD = 15.75*in2m;

% VESSEL (MAIN SECTIONS)
d0 = dump.OD;
d1 = 5.834*in2m;    % from 0 to box, excluding helicon
d2 = helicon.r*2;   % helicon
d3 = 24*in2m;       % ECH box
d4 = 19.25*in2m;    % Wide sections from box to end
d5 = 4.272*2*in2m;  % Narrow sections are just inside of coil
boxz1 = cmax(6);    % ECH box
boxz2 = cmin(7);
%           DUMP---------------------------DUMP,HELICON------------------------------HELICON,ECHBOX----------ECHBOX
vessel.z = [dump.z,dump.z,dump.z+0.2,dump.z+0.2,helicon.z1,helicon.z1,helicon.z2,helicon.z2,boxz1,boxz1,boxz2,boxz2];
vessel.r = [0     ,d0    ,d0         ,d1       ,d1        ,d2        ,d2        ,d1        ,d1   ,d3   ,d3   ,d5   ]/2;
for i=7:11
    vessel.z(end+1:end+4) = [cmax(i),cmax(i),cmin(i+1),cmin(i+1)];
    vessel.r(end+1:end+4) = [d5     ,d4     ,d4       ,d5       ]/2;
end
vessel.z(end+1) = 5;
vessel.r(end+1) = d5/2;

% OPTIONAL LIMITING STRUCTURES (skimmer, sleeve)
ves_add = [];
if add_skimmer
    ves_add.skimmer = skimmer;
end
if add_sleeve
    ves_add.sleeve = sleeve;
    kill_inds = find(vessel.z >= sleeve.z1 & vessel.z <= sleeve.z2);
    vessel.z(kill_inds) = [];
    vessel.r(kill_inds) = [];
end

% OVERLY COMPLICATED METHOD TO ADD COMPONENTS
if ~isempty(ves_add)
    names = fieldnames(ves_add);
    for i = 1:length(names)
        vessel.z = [vessel.z,ves_add.(names{i}).z];
        vessel.r = [vessel.r,ves_add.(names{i}).r.*ones(size(ves_add.(names{i}).z))];
    end
    [vessel.z,sortinds] = sort(vessel.z);
    vessel.r = vessel.r(sortinds);
    
    vessel.z2 = vessel.z(1);
    vessel.r2 = vessel.r(1);
    for i=2:length(vessel.z)-1
        if (vessel.r(i) ~= vessel.r(i-1)) && ~(vessel.z(i) == vessel.z(i-1))
            if vessel.r(i) > vessel.r(i-1)
                vessel.z2(end+1) = vessel.z(i-1);
                vessel.r2(end+1) = vessel.r(i);
            else
                vessel.z2(end+1) = vessel.z(i);
                vessel.r2(end+1) = vessel.r(i-1);
            end
        end
        vessel.z2(end+1) = vessel.z(i);
        vessel.r2(end+1) = vessel.r(i);
    end
    vessel.z2(end+1) = vessel.z(end);
    vessel.r2(end+1) = vessel.r(end);
    vessel.r = vessel.r2;
    vessel.z = vessel.z2;
    vessel = rmfield(vessel,'r2');
    vessel = rmfield(vessel,'z2');
end

% USED FOR TERMINATING FIELDLINES
vessel_clip_z = [vessel.z,5];
vessel_clip_r = [vessel.r,0];

% Combine into structure
geo.target = target;
geo.helicon = helicon;
geo.vessel =  vessel;
geo.skimmer = skimmer;
geo.coilcx.r = rcoil;
geo.coilcx.z = zcoil;
geo.cmax = cmax;
geo.cmin = cmin;
geo.vessel_clip_r = vessel_clip_r;
geo.vessel_clip_z = vessel_clip_z;
geo.sleeve = sleeve;
    

if plotit
    if newfig
        figure; hold on; box on;
    end
    % COILS
    for i = 1:size(zcoil,1)
        plot(zcoil(i,:),rcoil(i,:),'r')
        text(zcoil(i,1),mean(rcoil(i,1:2)),sprintf('%d',i),'color','r','fontsize',8,'fontweight','bold')
    end
    
    % TARGET
    plot(geo.target.z*[1,1],geo.target.r*[0,1],'k','linewidth',3)
    % HELICON
    plot([geo.helicon.z1,geo.helicon.z2],geo.helicon.r*[1,1],'k','linewidth',3)
    % VESSEL
    plot(geo.vessel.z,geo.vessel.r,'k','linewidth',1)    
    
%     plot(vessel_clip_z,vessel_clip_r,'c')
    
    xlabel('Z [m]','fontsize',14)
    ylabel('R [m]','fontsize',14)
    set(gca,'fontsize',14)
    axis([0,5,0,0.5])

end