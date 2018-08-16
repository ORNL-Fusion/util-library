function write_counter(ind,ind_max,text1)

% ic = abs(ind/ind_max*100-[0:10:100])<1;
% if any(ic) 
%     fprintf('%f\n',ind/ind_max*100)
% end
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