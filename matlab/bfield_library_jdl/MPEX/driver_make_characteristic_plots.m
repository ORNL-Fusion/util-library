clearvars;
shots = 7400 + [0,3:6,8,10,12:13,16,17,18]; mytitle = 'I_A = 6600 A, no skimmer';  x0_guess = -.5512;y0_guess = -2.533; force_guess = 1; scale_current_H = 2;

for i=1:length(shots)
    
    shot = shots(i);    
    fprintf('Working on shot %d, %d of %d\n',shot,i,length(shots))
    f = find_lcfs(shot,1);
    drawnow;
    make_characteristic_plot(shot);
    drawnow;
end