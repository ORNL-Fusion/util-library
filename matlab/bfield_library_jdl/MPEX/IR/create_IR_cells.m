function cells = create_IR_cells(IRdata,xshift,yshift)
dx = IRdata.dx;
dy = IRdata.dy;
nw = IRdata.nw;
nh = IRdata.nh;

x = linspace(-dx/2,dx/2,nw+1) - xshift;
y = linspace(-dy/2,dy/2,nh+1) - yshift;


icount = 1;
xcell = zeros(4,nh*nw);
ycell = zeros(4,nh*nw);
for i = 1:nw
    for j = 1:nh
        xcell(:,icount) = [x(i),x(i+1),x(i+1),x(i)];
        ycell(:,icount) = [y(j),y(j),y(j+1),y(j+1)];
        icount = icount + 1;
    end
end

xinterp1D = x(1:end-1) + diff(x)/2;
yinterp1D = y(1:end-1) + diff(y)/2;
[xmesh,ymesh] = meshgrid(xinterp1D,yinterp1D);

cells.x = x;
cells.y = y;
cells.xinterp1D = xinterp1D;
cells.yinterp1D = yinterp1D;
cells.xmesh = xmesh;
cells.ymesh = ymesh;
cells.ycell = ycell;
cells.xcell = xcell;