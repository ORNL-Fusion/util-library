function plot_coil_cross_section(rcoil,zcoil,newfig,current)


if nargin < 3
    newfig = 0;
end
if nargin < 4
    current = [];
end

if ~isequal(size(rcoil),size(zcoil))
    error('rcoil and zcoil must have the same size')
end
if ~isempty(current)
    if ~isvector(current)
        error('current must be a vector with one value per coil')
    end
    if length(current) ~= size(rcoil,1)
        error('length(current) must match the number of coil cross-sections')
    end
end

if newfig 
    figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14);
end

for i = 1:size(rcoil,1)
    plot(zcoil(i,:),rcoil(i,:),'k')
    if ~isempty(current)
        patch(zcoil(i,:),rcoil(i,:),current(i))
    end
end

if ~isempty(current)
    colorbar;
end
