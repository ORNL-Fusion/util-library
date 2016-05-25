function plot_sonnet(fname_in,plotit,newfig,col,lw,jxa_test,jxi_test)
% reads .sno file

% jxa_test = 54;
% jxi_test = 39;

if nargin < 2
    plotit = 1;
end
if nargin < 3
    newfig = 1;
end
if nargin < 4
    col = 'k';
end
if nargin < 5
    lw = 1;
end
if nargin < 6
    jxa_test = [];
end
if nargin < 7 
    jxi_test = [];
end
    
fid = fopen(fname_in,'r');
dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid); % thru R*Btor
dat = fgetl(fid); data = textscan(dat,'%s %s %d'); nx = cell2mat(data(3)); 
dat = fgetl(fid); data = textscan(dat,'%s %s %d'); ny = cell2mat(data(3));
dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid); dat = fgetl(fid);

ncell = (nx+2)*(ny+2);
fprintf('nx=%d, ny=%d, ncell=%d\n',nx,ny,ncell)
% Account for extra line in ddn grids
if isempty(dat) 
    dat = fgetl(fid);
end

icount = 1;
% for i = 1:ncell
for i = 1:ny+2
    for j = 1:nx+2
    dat = fgetl(fid); data = textscan(dat,'%s %d %s (%d %s %d %s %s %f,%f) %s %f,%f');
    crx(icount,1) = data{9}; cry(icount,1) = data{10};  crx(icount,2) = data{12};  cry(icount,2) = data{13};
    dat = fgetl(fid); data = textscan(dat,'%s %s %s %f %s %f,%f');
    crx(icount,5) = data{6}; cry(icount,5) = data{7};
    dat = fgetl(fid); data = textscan(dat,'%s %f,%f) %s %f,%f');
    crx(icount,3) = data{2}; cry(icount,3) = data{3};  crx(icount,4) = data{5};  cry(icount,4) = data{6};
    dat = fgetl(fid);
    crx2(i,j,1:4) = crx(icount,1:4);
    cry2(i,j,1:4) = cry(icount,1:4);
    icount = icount + 1;
    
    end
end

fclose(fid);

if newfig
    figure; hold on; box on;
end
if plotit
    for i = 1:ncell
        px = crx(i,[1,2,4,3,1]);
        py = cry(i,[1,2,4,3,1]);
        plot(px,py,[col,'-'],'linewidth',lw)
    end
end

if ~isempty([jxa_test,jxi_test])
    figure; hold on; box on;
    for i = 0:ny+1
        for j = 0:nx+1
            px = squeeze(crx2(i+1,j+1,[1,2,4,3,1]));
            py = squeeze(cry2(i+1,j+1,[1,2,4,3,1]));            
            if j == jxa_test
                plot(px,py,'r-','linewidth',3)
            elseif j == jxi_test 
                plot(px,py,'b-','linewidth',3)
            else
                plot(px,py,'k-','linewidth',2)
            end
        end
    end
end