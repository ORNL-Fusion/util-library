clearvars;
shots = 7400 + [0,3:6,8,10,12:13,16,17,18]; mytitle = 'I_A = 6600 A, no skimmer';
% shots = 7400 + [77,87,88,92:98,100,101,103]; mytitle = 'I_A = 3300 A, with skimmer';

i = 1

% for i=1:length(shots)    
    shot = shots(i);    
%     fprintf('Working on shot %d, %d of %d\n',shot,i,length(shots))
%     f = find_lcfs(shot,1);
    make_characteristic_plot(shot);
% end

