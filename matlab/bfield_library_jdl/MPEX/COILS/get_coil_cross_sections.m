function [rcoil,zcoil] = get_coil_cross_sections

fil = define_proto_coil_filaments;
rcoil = zeros(fil.ncoils,5);
zcoil = zeros(fil.ncoils,5);
for i = 1:fil.ncoils
    rcoil(i,:) = [fil.rr1(i),fil.rr2(i),fil.rr2(i),fil.rr1(i),fil.rr1(i)];
    zcoil(i,:) = [fil.z0(i),fil.z0(i),fil.z0(i)+fil.cl(i),fil.z0(i)+fil.cl(i),fil.z0(i)];
end

% figure; hold on; box on;
% for i = 1:ncoils
%     plot(zcoil(i,:),rcoil(i,:),'r')
% end