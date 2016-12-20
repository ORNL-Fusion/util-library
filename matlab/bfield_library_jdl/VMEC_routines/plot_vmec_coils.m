function plot_vmec_coils(coil,nfp)

if nargin < 2
    nfp = 1;
end
nc = length(coil.coilpos);
figure; hold on; box on;
for i = 1:nc
    coilx = coil.coilpos{i}(:,1);
    coily = coil.coilpos{i}(:,2);
    coilz = coil.coilpos{i}(:,3);
    if nfp > 1
        coilp = atan2(coily,coilx);
        coilr = sqrt(coilx.^2 + coily.^2);
        [coilr,coilp,coilz] = symmetrize_vec(coilr,coilp,coilz,nfp,0);
        coilx = coilr.*cos(coilp);
        coily = coilr.*sin(coilp);
    end
    plot3(coilx,coily,coilz)
    
end
