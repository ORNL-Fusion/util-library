function flare_plot_poincare(fname,newfig,sortit,plotstr)
  if nargin < 1
    error('Must specify fname')    
  end
  if nargin < 2
    newfig = 1;
  end
  if nargin < 3
    sortit = 1;
  end
  if nargin < 4
    plotstr = '.-';
  end
  
%  fname = fullfile(run_path,'poincare_vac_0.dat');
  
  if newfig 
    figure; hold on; box on;
  end

  fprintf('Reading file: %s\n',fname)

  %% Check for comment lines
  iskip = 0;
  fid = fopen(fname,'r');
  dat = fgetl(fid);
  while strcmp(dat(1),'#')
    iskip = iskip + 1;
    dat = fgetl(fid);
  end  
  fclose(fid);
  iskip;

  d = dlmread(fname,'',iskip,0);
  R     = d(:,1);
  Z     = d(:,2);
  if size(d,2) > 2
    theta = d(:,3);
  end
  if size(d,2) > 3
    psi   = d(:,4);
  end
  if size(d,2) > 4
    ind   = d(:,5);
  end
  
  if ~sortit
    plot(R,Z,plotstr);
    return
  end

  fstarts = find(diff(ind) <= 0) + 1;
  istarts = [1;fstarts];
  iends   = [fstarts - 1;length(R)];
  N_surfs = length(istarts);
  for i = 1:N_surfs
    i1 = istarts(i);
    i2 = iends(i);
%    fprintf('inds : %d %d\n',ind(i1),ind(i2))
    [~,isort] = sort(theta(i1:i2));
    Rtmp = R(i1:i2);
    Ztmp = Z(i1:i2);
    Rsort{i} = [Rtmp(isort);Rtmp(isort(1))];
    Zsort{i} = [Ztmp(isort);Ztmp(isort(1))];
    plot(Rsort{i},Zsort{i},plotstr,'markersize',10);
  end

  

  
