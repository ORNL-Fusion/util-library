function ts2d = read_2dts_dat_file(fname)
% CHANNEL INDICES ARE OFFSET BY 1!!!!
% fname = 'C:\Work\DIII-D\164723\2DTS_Stuff\dts_fitted1d_osp_164722_2600to5000msEFIT01_remapto141905at3500msEFIT01_fuelD_ioncharge1_gamma7p539_see0p000.dat';


fid = fopen(fname);

fgetl(fid);
fgetl(fid);
fgetl(fid);

ichan_last = -1;
icount = 0;
while ~feof(fid)
% d=fgetl(fid); 
% numc = find(d == ',');

d1 = fscanf(fid,'%d %f',2);  % shot, laser time
d2 = fscanf(fid,'%s',1); % Chan
d3 = fscanf(fid,'%f',2); % R,dts [m] Z,dts [m]
d4 = fscanf(fid,'%s',1); % EFIT RUNID
d5 = fscanf(fid,'%f',2); % EFIT time, PSIn
d6 = fscanf(fid,'%s',1); % Target region
d7 = fscanf(fid,'%f',6); % Rx, Zx, Rvsout, Zvsout, rvsin, zvsin
d8 = fscanf(fid,'%f',8); % Lpol xpt to OSP, to ISP, DTS to OSP, DTS to ISP, 4 more Ls
d9 = fscanf(fid,'%f',7); % Rremap, Zremap, Rmaptarg, Zmaptarg, R-Rsep map, Z-Zsep map, ROMP
d10= fscanf(fid,'%f',5); %Te, Teerr, err ratio, Te lowess wid, TEfit
d11= fscanf(fid,'%f',5); % same for ne [10^20]
d12= fscanf(fid,'%f',5); % same for pe [Pa]
d13= fscanf(fid,'%f',4); % Ti/Te ratio, Ti, Ti err, err rat
d15= fscanf(fid,'%s',1); % Fuel ion species
d16= fscanf(fid,'%f\n',34); % stuff

ichan = str2num(d2);

if length(ichan) ~= 1
    continue;
else  
    if ichan ~= ichan_last
        icount = 0;
    end
    ichan_last = ichan;
    icount = icount + 1;
    ts2d.chan{ichan+1}.Te(icount) = d10(1);
    ts2d.chan{ichan+1}.dTe(icount) = d10(2);
    ts2d.chan{ichan+1}.ne(icount) = d11(1);
    ts2d.chan{ichan+1}.dne(icount) = d11(2);    
    ts2d.chan{ichan+1}.Romp(icount) = d9(7);
    ts2d.chan{ichan+1}.psiN(icount) = d5(2);    
    ts2d.chan{ichan+1}.Tefit(icount) = d10(5);
    ts2d.chan{ichan+1}.nefit(icount) = d11(5);
end

end

% fclose(fid);