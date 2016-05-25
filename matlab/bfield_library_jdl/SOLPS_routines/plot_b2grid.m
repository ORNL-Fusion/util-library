function plot_b2grid(fname_in,plotit,newfig,col,lw)
% This will read the .geo file
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

data = dlmread(fname_in);

nr = data(1,1) + 2; 
np = data(1,2) + 2;

icount = 2;
for i = 1:nr
    for j = 1:np
        grid(i,j,:) = data(icount,3:14);
        xgrid(i,j,:) = data(icount,3:2:12);   % First is cell center, then 4 corners
        ygrid(i,j,:) = data(icount,4:2:12);
        icount = icount + 1;
    end
end

if newfig
    figure; hold on; box on;
end
if plotit
    for i = 1:nr
        for j = 1:np
            px = squeeze(xgrid(i,j,[2,3,5,4,2]));
            py = squeeze(ygrid(i,j,[2,3,5,4,2]));
            plot(px,py,[col,'-'],'linewidth',lw)
        end
    end
end


% figure; hold on; box on;
% plot(data(2:end,3),data(2:end,4),'o')
