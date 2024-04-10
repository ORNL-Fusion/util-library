function f = get_resScale(nxPlot,nyPlot,gap,botBuffer,topBuffer,leftBuffer,rightBuffer)
%f = get_resScale(nxPlot,nyPlot,gap,botBuffer,topBuffer,leftBuffer,rightBuffer)

resRaw = get(0,'ScreenSize');
xScale = resRaw(3)/3440;
yScale = resRaw(4)/1440;
f.resScale = [xScale yScale 1 1];


if nargin == 0
    return
end
% 
% if nxPlot ~= 1
%     error('Update this for nx > 1')
% end

if nargin < 3
    gap = 0.04;
end
if nargin < 4
    botBuffer = 0.15;
end
if nargin < 5
    topBuffer = 0.05;
end
if nargin < 6
    leftBuffer = 0.19;
end
if nargin < 7
    rightBuffer = 0.05;
end



% dx = 1 - leftBuffer - rightBuffer;
yMin = botBuffer;
yMax = 1 - topBuffer;
xMin = leftBuffer;
xMax = 1 - rightBuffer;

dy = (yMax - yMin - gap*(nyPlot - 1))/nyPlot;
dx = (xMax - xMin - gap*(nxPlot - 1))/nxPlot;


for iy = 1:nyPlot
    for ix = 1:nxPlot
        % iPlot = iy + (ix-1)*nyPlot;
        % iPlot = ix + (iy-1)*nxPlot;
        f.p(ix,iy,:) = [(dx+gap)*(ix-1) + xMin,(dy+gap)*(nyPlot-iy) + yMin,dx,dy];
    end
end