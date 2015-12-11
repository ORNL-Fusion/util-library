function plot_c2_grid(grid)

cf = colorflipper(grid.ndomain,'jet');
figure; hold on; box on;
% set up legend
for idomain = 1:grid.ndomain
   plot(NaN,'-','color',cf(idomain,:))
   leg_str{idomain} = num2str(idomain);
end
for idomain = 1:grid.ndomain
    plot(grid.x2d{idomain},grid.y2d{idomain},'-','color',cf(idomain,:))
    plot(grid.x2d{idomain}.',grid.y2d{idomain}.','-','color',cf(idomain,:))
end
legend(leg_str);