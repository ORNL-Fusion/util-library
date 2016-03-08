function data = load_IR_file(fname,plotit)
% Note that image is fliped up-down
if nargin < 2
    plotit = 0;
end

px_per_cm = 12.146;
d = load(fname);

% plot raw data
if plotit
    figure; hold on; box on;
    h = pcolor(d.Frame);
    colorbar;
    set(h,'edgecolor','none')
    title('Raw data','fontsize',14)
    xlabel('X [px]','fontsize',14)
    ylabel('Y [px]','fontsize',14)
    set(gca,'fontsize',14)
end

% Convert to cm and flip
IRdata2D = flipud(d.Frame);
nw = size(d.Frame,2);
nh = size(d.Frame,1);
dx = nw/px_per_cm;
dy = nh/px_per_cm;
IRdata1D = reshape(IRdata2D,1,nw*nh);

data.dx = dx;
data.dy = dy;
data.nw = nw;
data.nh = nh;
data.IRdata1D = IRdata1D;
data.IRdata2D = IRdata2D;
