function elm=check_elm_timing(elm_file,twin)
% clearvars;

myelm = load(elm_file);
% load(fs_file);

% twin = [3000,4900];
elm_norm_thresh = 0.4;
elm_dt_thresh = 10;

igood = find(myelm.elm.telmpeak >= twin(1) & myelm.elm.telmpeak <= twin(2));
telm = myelm.elm.telmpeak(igood);
elmpeak = myelm.elm.elmpeak(igood);

igood_fs = find(myelm.elm.t >= twin(1) & myelm.elm.t <= twin(2));
tfs = myelm.elm.t(igood_fs);
fs = myelm.elm.d(igood_fs);


emean = mean(elmpeak);
elmpeak_norm = elmpeak/max(elmpeak);

ielm_keep = find(elmpeak_norm >= elm_norm_thresh);


tkeep = telm(ielm_keep);



% dt filter
dtkeep = diff(tkeep);
ikeep_dt = find(dtkeep >= elm_dt_thresh);
ielm_keep = ielm_keep(ikeep_dt);

% Final results
tkeep = tkeep(ikeep_dt);
nkeep = length(tkeep);

fs_at_elm = interp1(tfs,fs,tkeep);

figure; hold on;
plot(tfs,fs./max(fs),'k')
plot(telm,elmpeak_norm,'r.')
% plot(telm(ielm_keep),elmpeak_norm(ielm_keep),'ro')
plot(tkeep,fs_at_elm./max(fs),'ro')
% plot(telm,emean*ones(size(telm)),'k')
ylabel('Norm. filterscope')
xlabel('Time [ms]')
title('ELM timing')
legend('fs','elm marks','selected')

elm.telm = tkeep;
elm.nelm = nkeep;
