function plot_c2_grid(grid)

cf = colorflipper(grid.ndomain,'jet');
figure; hold on; box on;
for idomain = 1:grid.ndomain
    plot(grid.x2d{idomain},grid.y2d{idomain},'-','color',cf(idomain,:))
    plot(grid.x2d{idomain}.',grid.y2d{idomain}.','-','color',cf(idomain,:))
end