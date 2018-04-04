function plot_w7x_IR_dynamic_clim(image,min_guess,max_guess)
    if nargin < 2
        min_guess = 0;
    end
    if nargin < 3
        max_guess = max(image(:));
    end
    % clearvars;
    
    % load('C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR\shot27_irdata_cam10_t2_3_image1.mat')
    
%     min_guess = 350;
%     max_guess = 500;
    
    figure; hold on; box on;
    pp = pcolor(image);
    set(pp,'edgecolor','none')
        set(gca,'clim',[min_guess,max_guess])
  
    cb = colorbar;
    axis tight;
    
    Min = max(0,min(image(:)));
    Max = max(max(image(:)),max_guess);
    Pos=get(cb,'position');
    
    uic_min=uicontrol('style','slider','units','normalized',...
        'position',[Pos(1)+0.11 Pos(2) 0.04 0.4],...
        'min',Min,'max',Max,'Value',min_guess,...
        'callback',@Smin);
    uic_max=uicontrol('style','slider','units','normalized',...
        'position',[Pos(1)+0.11 0.5 0.04 0.4],...
        'min',Min,'max',Max,'Value',max_guess,...
        'callback',@Smax);
    
    function Smin(source,~)
        cc = get(gca,'clim');
        set(gca,'clim',[source.Value,cc(2)]);
        drawnow;
    end
    function Smax(source,~)
        cc = get(gca,'clim');
        set(gca,'clim',[cc(1),source.Value]);
        drawnow;
    end
    
end