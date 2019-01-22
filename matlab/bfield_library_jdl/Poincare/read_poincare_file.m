function poinc = read_poincare_file(run_path,fname,fname2,g)



for ifile = 1:2
    if ifile == 1
        fname_tmp = fname;
    else
        fname_tmp = fname2;
    end
    fid = fopen([run_path,fname_tmp]);
    if fid == -1
        error('ERROR: cannot open %s\n',[run_path,fname_tmp])
    end
    dat = sscanf(fgets(fid),'%f %i %i',3);
    phi_plot_deg = dat(1);
    numlines = dat(2); poinc.numlines = numlines;
    num_pts = dat(3);
    
    rline{ifile} = NaN(numlines,num_pts);
    zline{ifile} = NaN(numlines,num_pts);
    psiNline{ifile} = NaN(numlines,num_pts);
    iline{ifile} = NaN(numlines,1);
    for i = 1:numlines
        iline{ifile}(i) = fscanf(fid,'%i',1);
        try
        rline{ifile}(i,:) = fscanf(fid,'%f',num_pts);
        catch
            a=1;
        end
        zline{ifile}(i,:) = fscanf(fid,'%f',num_pts);
        if ~isempty(g)
            psiNline{ifile}(i,:) = calc_psiN(g,rline{ifile}(i,:),zline{ifile}(i,:),1);
        end
    end
    fclose all;
    
    %     itmp = isnan(psiNline{ifile});
    
    if ~isempty(g)
        rline{ifile}(isnan(psiNline{ifile})) = NaN;
        zline{ifile}(isnan(psiNline{ifile})) =NaN;
        psiNline{ifile}(isnan(psiNline{ifile})) = NaN;
    end
end

poinc.rline = rline;
poinc.zline = zline;
poinc.psiNline = psiNline;