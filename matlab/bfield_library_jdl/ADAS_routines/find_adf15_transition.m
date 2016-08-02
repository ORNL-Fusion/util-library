function [i_ex,i_rc,i_cx] = find_adf15_transition(adf15,wl_nm,plotit)
if nargin < 3
    plotit = 0;
end
% fname = 'C:\Work\ADAS\pec12#h_pju#h0.dat';
% adf15 = read_adas_adf15_file(fname);
% wl_nm = 656;

dl_ang_tol = 10;
[dl_ang,ind_min] = min(abs(adf15.wlng-wl_nm*10));
wl_found_ang = adf15.wlng(ind_min);

both_inds = find(abs(adf15.wlng - wl_found_ang) < 1e-6);
if length(both_inds) < 2 || length(both_inds) > 3
    error('Only 2 or 3 reactions are understood')
end    

i_cx = [];
for i = 1:length(both_inds)
    switch adf15.type{both_inds(i)}
        case 'EXCIT'
            i_ex = both_inds(i);
        case 'RECOM'
            i_rc = both_inds(i);
        case 'CHEXC'
            i_cx = both_inds(i);            
        otherwise
            error('Did not recognize type')
    end
end

fprintf('-----------------------------------------------\n')
fprintf('Wanted wl = %8.3f [nm], found %8.3f [nm]\n',wl_nm,wl_found_ang/10)
if isempty(i_cx)
    fprintf('i_ex = %d, i_rc = %d\n',i_ex,i_rc)
else
    fprintf('i_ex = %d, i_rc = %d, i_cx = %d\n',i_ex,i_rc,i_cx)
end
fprintf('-----------------------------------------------\n')
if abs(dl_ang) > dl_ang_tol
    error('wl deviation greater than tol')
end

if plotit
    figure; hold on; box on;
    plot(adf15.te{i_ex},adf15.pec{i_ex},'-')
    plot(adf15.te{i_rc},adf15.pec{i_rc},'--')
    if ~isempty(i_cx)
        plot(adf15.te{i_cx},adf15.pec{i_cx},'--')
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
