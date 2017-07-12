function xdr = read_xdr_dump_file(out_path,fname)
% clearvars;

% out_path = 'C:\Work\Stellarator\ALL_W7X_WORK\xdr_dump_read\OUTPUT\';
% fname = 'field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.out';

% check for mat file, create it if it does not exist
fn_mat = [out_path,fname(1:end-3),'mat'];
if exist(fn_mat,'file') ~= 2
    fn_out = [out_path,fname];
    disp('Loading data from .out file and writing .mat file');
    fid = fopen(fn_out,'r');
    i1 = fscanf(fid,'%i',1);
    i2 = fscanf(fid,'%i',1);
    rnull = fscanf(fid,'%f',1);
    ronull = fscanf(fid,'%f',1);
    eta = fscanf(fid,'%f',1);
    fnull = fscanf(fid,'%f',1);
    nperio = fscanf(fid,'%i',1);
    ialfa = fscanf(fid,'%i',1);
    knull = fscanf(fid,'%i',1);
    iganz = fscanf(fid,'%i',1);
    ampby0 = fscanf(fid,'%f',1);
    bz0 = fscanf(fid,'%f',1);
    bfak = fscanf(fid,'%f',1);
    k2 = fscanf(fid,'%i',1);
    iald21 = fscanf(fid,'%i',1);
    
    brg = zeros(iald21,k2,k2);
    bfg = zeros(iald21,k2,k2);
    bzg = zeros(iald21,k2,k2);
    
    disp('Reading Bgrid info (probably going to take forever)')
    for iz=1:k2
        for ir=1:k2
            for i=1:iald21
                brg(i,ir,iz) = fscanf(fid,'%f',1);
            end
        end
    end
    disp('Done with brg')
    for iz=1:k2
        for ir=1:k2
            for i=1:iald21
                bfg(i,ir,iz) = fscanf(fid,'%f',1);
            end
        end
    end
    disp('Done with bfg')
    for iz=1:k2
        for ir=1:k2
            for i=1:iald21
                bzg(i,ir,iz) = fscanf(fid,'%f',1);
            end
        end
    end    
    disp('Done with bzg')
    
    fclose(fid);
    
    disp('Saving .mat file')
    save(fn_mat);
    
else
    disp('Loading data from .mat file');
    load(fn_mat);
end

xdr.fname = fname;
xdr.i1 = i1;
xdr.i2 = i2;
xdr.rnull = rnull;
xdr.ronull = ronull;
xdr.eta = eta;
xdr.fnull = fnull;
xdr.nperio = nperio;
xdr.ialfa = ialfa;
xdr.knull = knull;
xdr.iganz = iganz;
xdr.ampby0 = ampby0;
xdr.bz0 = bz0;
xdr.bfak = bfak;
xdr.k2 = k2;
xdr.iald21 = iald21;
xdr.brg = brg;
xdr.bfg = bfg;
xdr.bzg = bzg;


if iganz ~= 0
    error('I think this means no stellarator symmetry!')
end
if iganz ~= 0 || ampby0 ~= 0 || bz0 ~= 0 || fnull ~= 1
    error('old variables I have assumed have these values')
end