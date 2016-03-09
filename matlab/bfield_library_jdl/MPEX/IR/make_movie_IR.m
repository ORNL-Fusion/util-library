clearvars;
% shots = 7400 + [0,3:6,8,10,12:13,16,17,18]; mytitle = 'I_A = 6368 A, no skimmer';  x0_guess = -.5512;y0_guess = -2.533; force_guess = 1;
% shots = 7400 + [3,8,12,18]; mytitle = 'I_A = 3300 A, no skimmer';
shots = 7400 + [77,87,88,92:98,100,101,103]; mytitle = 'I_A = 3300 A, with skimmer';x0_guess = -.5877;y0_guess = -2.8914; force_guess = 1;
ns = length(shots);
% x0_guess = -.4;y0_guess = -3.2;



scale_H = 1;

plotitIR = 1;
for i = 1:ns
    shot = shots(i);
    [helicon_current(i),current_A,current_B,config,skimmer] = get_Proto_current(shot);
    [rr{i},dd{i},radius(i),angle(i),x0_final(i),y0_final(i)] = fit_IR_data(shot,plotitIR,x0_guess,y0_guess,force_guess);
    set(gca,'clim',[0,14])
    axis([-4,4,-4,4])
    drawnow;
    F(i) = getframe;
    aa = 1
end


fig = figure;
movie(fig,F)
