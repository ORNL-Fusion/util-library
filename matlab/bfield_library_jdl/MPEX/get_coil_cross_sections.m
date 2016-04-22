function [rcoil,zcoil] = get_coil_cross_sections

[nturns,nlayers,rr1,rr2,cl,z0] = define_proto_coil_filaments;
ncoils = length(nturns);
for i = 1:ncoils
    rcoil(i,:) = [rr1(i),rr2(i),rr2(i),rr1(i),rr1(i)];
    zcoil(i,:) = [z0(i),z0(i),z0(i)+cl(i),z0(i)+cl(1),z0(i)];
end

% figure; hold on; box on;
% for i = 1:ncoils
%     plot(zcoil(i,:),rcoil(i,:),'r')
% end