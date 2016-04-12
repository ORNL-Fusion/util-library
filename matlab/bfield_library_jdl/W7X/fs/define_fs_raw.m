function fs = define_fs_raw(plotit,newfig)
if nargin < 1
    plotit = 0;
    newfig = 0;
end
if nargin == 1
    newfig = 1;
end

% AEA30
p1 = [-6214.195,4514.877,-1057.225]/1000;
p2 = [-4614.658,3352.745,84.275]/1000;
fs.AEA30.p1 = p1;
fs.AEA30.p2 = p2;

% AEK41
p1 = [-5098.138,-6557.635,369.232]/1000;
p2 = [-5034.268,-6480.69,369.2]/1000;
fs.AEK41.p1 = p1;
fs.AEK41.p2 = p2;


% % AEK51 turned by 180
% p1 = [4607.795,-6790.411,369.184]/1000;
% p2 = []/1000;

% AEK51
p1 = [4390.861,-7046.029,369.184]/1000;
p2 = [4337.388,-6961.527,369.184]/1000;
fs.AEK51.p1 = p1;
fs.AEK51.p2 = p2;

% AEL41
p1 = [-2263.798,-2072.288,-160.018]/1000;
p2 = [-2308.574,-2160.724,-146.815]/1000;
fs.AEL41.p1 = p1;
fs.AEL41.p2 = p2;

% AEL51
p1 = [1340.786,-2860.601,-145.792]/1000;
p2 = [1411.057,-2930.514,-132.589]/1000;
fs.AEL51.p1 = p1;
fs.AEL51.p2 = p2;

if newfig
    figure; hold on; box on;
end
if plotit
    names = fieldnames(fs);
    nfs = length(names);    
    for i = 1:nfs
        ifs = getfield(fs,names{i});
        plot3([ifs.p1(1),ifs.p2(1)],[ifs.p1(2),ifs.p2(2)],[ifs.p1(3),ifs.p2(3)])
    end
end
    
