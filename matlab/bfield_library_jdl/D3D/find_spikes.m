function spikes = find_spikes(t,d,thresh,dt_thresh,MEDFILT)

if MEDFILT
d=d-medflt1d_jl(d,floor(length(t)*0.1),'zero');
end

med = median(d);
sig_norm = (d - med)./max(d-med);


% Find peaks
right_neighbor = [0;sig_norm(1:end-1)];
left_neighbor = [sig_norm(2:end);0];
p_ind = find((sig_norm > right_neighbor) & (sig_norm > left_neighbor) & (sig_norm > thresh));

% eliminate multiples
dt_ind = diff(t(p_ind));
max_loop = 10;
iloop = 0;
while any(dt_ind < dt_thresh)

    i = 1;
    while i < length(dt_ind) + 1
        if dt_ind(i) < dt_thresh %&& dt_ind(i+1) < dt_thresh
            if sig_norm(p_ind(i)) > sig_norm(p_ind(i+1))
                p_ind(i+1) = [];
                i = i - 1;
            else
                p_ind(i) = [];
                i = i - 1;
            end
        end
        dt_ind = diff(t(p_ind));
        i = i+1;
    end
    
    iloop = iloop + 1;
    if iloop > max_loop
        disp(dt_ind)
        error('too many iterations')
    end    
end

figure; hold on; box on;
subplot(2,1,1); hold on; box on;
plot(t,d)
yline(med,'k')
subplot(2,1,2); hold on; box on;
ylabel('Heat flux [MW/m^2]')
xlabel('Time [ms]')
plot(t,sig_norm)
plot(t(p_ind),sig_norm(p_ind),'ro')
yline(thresh,'k')
ylabel('Normalized signal')
xlabel('Time [ms]')
title('ELM timing (IR)')

spikes.time_inds = p_ind;
spikes.times = t(p_ind);
spikes.val = d(p_ind);
spikes.val_norm = sig_norm(p_ind);