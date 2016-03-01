function [geo] = get_Proto_geometry(plotit,newfig)

if nargin < 1
    plotit = 0;
end
if nargin < 2
    newfig = 1;
end

add_skimmer = 1;

in2m = 0.0254;

[rcoil,zcoil] = get_coil_cross_sections;
cmax = max(zcoil,[],2);  % max and min Z of the coils
cmin = min(zcoil,[],2);


target.z = cmax(7) + (cmin(8) - cmax(7))/2 + 1*in2m;  % midway between 7-8, 1" towards 8
target.r = 3.75*in2m/2;  % diameter is 3.75 inches

helicon_ID = 5.1*in2m;                            % diameter is 5.1 inches (5.1 from Meitner)
helicon_L  = 11.8*in2m;                           % L is 11.8" from Meitner
helicon_Z_cen = cmax(3) + (cmin(4) - cmax(3))/2;  % Centered between coils 3 and 4
helicon.z1 = helicon_Z_cen - helicon_L/2;         
helicon.z2 = helicon_Z_cen + helicon_L/2;
helicon.r = helicon_ID/2; 
skimmer.z1 = (cmin(5) + cmax(5))/2 + 0.0702 - 0.5e-2;  
skimmer.z2 = (cmin(5) + cmax(5))/2 + 0.0702 + 0.5e-2;  
skimmer.ID = 0.058;

d1 = 5.834*in2m;  % from 0 to box, excluding helicon
d2 = helicon_ID;  % helicon
d3 = 24*in2m;     % box
d4 = 19.25*in2m;  % wide sections from box to end
d5 = 4.272*2*in2m;  % Narrow sections are just inside of coil
d6 = skimmer.ID;       % skimmer
boxz1 = cmax(6);
boxz2 = cmin(7);

if add_skimmer
    vessel.r = [d1,d1        ,d2        ,d2        ,d1        ,d1            ,d6            ,d6            ,d1            ,d1   ,d3   ,d3   ,d5   ,d5     ,d4     ,d4     ,d5     ,d5     ,d4     ,d4     ,d5     ,d5     ,d4     ,d4      ,d5      ,d5      ,d4      ,d4      ,d5      ,d5      ,d4      ,d4      ,d5      ,d5]/2;
    vessel.z = [0 ,helicon.z1,helicon.z1,helicon.z2,helicon.z2,skimmer.z1,skimmer.z1,skimmer.z2,skimmer.z2,boxz1,boxz1,boxz2,boxz2,cmax(7),cmax(7),cmin(8),cmin(8),cmax(8),cmax(8),cmin(9),cmin(9),cmax(9),cmax(9),cmin(10),cmin(10),cmax(10),cmax(10),cmin(11),cmin(11),cmax(11),cmax(11),cmin(12),cmin(12),5];
else
    vessel.r = [d1,d1        ,d2        ,d2        ,d1        ,d1   ,d3   ,d3   ,d5   ,d5     ,d4     ,d4     ,d5     ,d5     ,d4     ,d4     ,d5     ,d5     ,d4     ,d4      ,d5      ,d5      ,d4      ,d4      ,d5      ,d5      ,d4      ,d4      ,d5      ,d5]/2;
    vessel.z = [0 ,helicon.z1,helicon.z1,helicon.z2,helicon.z2,boxz1,boxz1,boxz2,boxz2,cmax(7),cmax(7),cmin(8),cmin(8),cmax(8),cmax(8),cmin(9),cmin(9),cmax(9),cmax(9),cmin(10),cmin(10),cmax(10),cmax(10),cmin(11),cmin(11),cmax(11),cmax(11),cmin(12),cmin(12),5];    
end
    

geo.target = target;
geo.helicon = helicon;
geo.vessel =  vessel;
geo.skimmer = skimmer;

if plotit
    if newfig
    figure; hold on; box on;
    end
    % COILS
    for i = 1:size(zcoil,1)
        plot(zcoil(i,:),rcoil(i,:),'r')
    end
    
    % TARGET
    plot(geo.target.z*[1,1],geo.target.r*[0,1],'k','linewidth',3)
    % HELICON
    plot([geo.helicon.z1,geo.helicon.z2],geo.helicon.r*[1,1],'k','linewidth',3)
    % VESSEL
    plot(geo.vessel.z,geo.vessel.r,'k','linewidth',1)
    
    xlabel('Z [m]','fontsize',14)
    ylabel('R [m]','fontsize',14)
    set(gca,'fontsize',14)
    axis([0,5,0,0.5])
    
end