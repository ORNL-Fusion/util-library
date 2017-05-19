function flare_plot_boundary(run_path,newfig)
  if nargin < 1
    run_path = '.';
  end
  if nargin < 2
    newfig = 1;
  end
  
  files = dir(fullfile(run_path,'boundary_quad_*.plt'));
  
  nf = length(files);
  fprintf('Found %d boundary_quad files\n',nf)
  
  if newfig 
    figure; hold on; box on;
  end
  
  for i = 1:nf
    d = dlmread(fullfile(run_path,files(i).name),'',1,0);
    plot(d(:,1),d(:,2),'k.-');
  end
  

