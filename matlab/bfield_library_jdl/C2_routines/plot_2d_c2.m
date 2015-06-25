function plot_2d_c2(grid,data_plot,mytitle)
if nargin < 3
    mytitle = ' ';
end

icount = 1;
for idomain = 1:grid.ndomain
    for i = 1:grid.ny(idomain)-1
        for j = 1:grid.nx(idomain)-1
             xc(:,icount) = [grid.x2d{idomain}(i,j),grid.x2d{idomain}(i,j+1),grid.x2d{idomain}(i+1,j+1),grid.x2d{idomain}(i+1,j)];
             yc(:,icount) = [grid.y2d{idomain}(i,j),grid.y2d{idomain}(i,j+1),grid.y2d{idomain}(i+1,j+1),grid.y2d{idomain}(i+1,j)];
             d(icount) = data_plot{idomain}(i,j);
             icount = icount + 1;  
        end
    end
end

figure; hold on; box on; 
patch(xc(:,~isnan(d)),yc(:,~isnan(d)),d(~isnan(d)),'edgecolor','none')
fprintf('Max value in this plot: %f\n',max(d(~isnan(d))))

axis tight;
colorbar('fontsize',14)
xlabel('R [m]','fontsize',14)
ylabel('Z [m]','fontsize',14)
set(gca,'fontsize',14)
title(mytitle);

colorbar;