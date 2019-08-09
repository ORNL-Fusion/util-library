function write_counter(ind,ind_max,text1)
% Creates a labeled counter that deletes the previous entry. 
% Ex: 
% for i = 1:10
%     write_counter(i,10,"Number of surfaces to compute:")
%     pause(0.25);
% end
% 
% Will print:
% Number of surfaces to compute: 10; working on: 10
% Where the last number will increment from 1 to 10 in-line

if ind > 1
    for iloop = 0:log10(ind-1)
        fprintf('\b');
    end
    fprintf('%d',ind);
    if ind == ind_max
        fprintf('\n')
    end
else
    fprintf('  %s %d; working on: %d',text1,ind_max,ind);
end