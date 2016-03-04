function fout = clip_fl_at_vessel(f,vessel_r,vessel_z)
if nargin == 2    % then vessel_r should be geo structure
    vessel_z = vessel_r.vessel_clip_z;
    vessel_r = vessel_r.vessel_clip_r;
end
for i=1:size(f.r,2)
    isin2 = inpolygon(f.r(:,i),f.z(:,i),vessel_r,vessel_z);
    is2 = find(isin2 == 0,1,'first')-1;
    if isempty(is2)
        is2 = length(isin2);
    end
    fout.r(:,i) = [f.r(1:is2,i);NaN(length(is2+1:size(f.r,1)),1)];
    fout.z(:,i) = [f.z(1:is2,i);NaN(length(is2+1:size(f.r,1)),1)];
end
