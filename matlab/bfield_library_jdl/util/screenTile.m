function Tile = screenTile(nCol,nRow)


Tile.nCol = nCol;
Tile.nRow = nRow;

screenPx = get(0,'screensize');
pxHoz = screenPx(3);
pxVer = screenPx(4);
hPxBuf = floor(max(pxHoz*0.003,10));
vPxBuf = floor(max(pxVer*0.07,100));



Tile.dHoz = floor((pxHoz - 2*hPxBuf - hPxBuf*(Tile.nCol-1))/Tile.nCol);
Tile.dVer = floor((pxVer - 2*vPxBuf - vPxBuf*(Tile.nRow-1))/Tile.nRow);

Tile.hPosStart = hPxBuf + (Tile.dHoz+hPxBuf)*[0:Tile.nCol-1];
Tile.vPosStart = vPxBuf + (Tile.dVer+vPxBuf)*[0:Tile.nRow-1];



% for i=1:Tile.nCol
%     for j = 1:Tile.nRow
%          figure('Position',[Tile.hPosStart(i),Tile.vPosStart(j),Tile.dHoz,Tile.dVer])
%     end
% end
    
