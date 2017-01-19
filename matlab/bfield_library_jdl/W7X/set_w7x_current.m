function coil_or_mgrid = set_w7x_current(coil_or_mgrid,taper)
           %1  2  3  4  5  6  7   8   9       10       11       12       13       14
% taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2,Itrim_A1,Itrim_A2,Itrim_A3,Itrim_A4,Itrim_B1];
% taper should be the current per winding, as it is multiplied by the number of windings in each coil below
taper_tmp = taper;
taper = zeros(1,14);
taper(1:length(taper_tmp)) = taper_tmp;


% Figure out if mgrid or coil
iscoil = 0;
ismgrid = 0;
if isfield(coil_or_mgrid,'cname')
    iscoil = 1;
    nc = length(coil_or_mgrid.coilpos);
elseif isfield(coil_or_mgrid,'coil_group')
    ismgrid = 1;
    nc = length(coil_or_mgrid.raw_coil_cur);
else
    error('Did not recognize what coil_or_mgrid is')
end



current = [];
for i = 1:nc
    if iscoil
        cname = strtrim(lower(char(coil_or_mgrid.cname{i})));
    end
    if ismgrid
        cname = strtrim(lower(coil_or_mgrid.coil_group(:,i).'));
    end
    switch cname
        case 'coil_1'
            tval = taper(1);
        case 'coil_2'
            tval = taper(2);
        case 'coil_3'
            tval = taper(3);
        case 'coil_4'
            tval = taper(4);
        case 'coil_5'
            tval = taper(5);
        case 'coil_a'
            tval = taper(6);
        case 'coil_b'
            tval = taper(7);
        case 'sweep_coil_1'
            tval = taper(8);
        case 'sweep_coil_2'
            tval = taper(9);
        case 'trim_a1'
            tval = taper(10);
        case 'trim_a2'
            tval = taper(11);
        case 'trim_a3'
            tval = taper(12);
        case 'trim_a4'
            tval = taper(13);
        case 'trim_b1'
            tval = taper(14);            
        otherwise
            error('Did not recognize coil type in set_w7x_current %s\n',cname)
    end    
    if iscoil
        numwind = coil_or_mgrid.numwind(i);
        coilcur{i} = tval*ones(size(coil_or_mgrid.coilcur{i}))*numwind;
        coilcur{i}(end) = 0;
        current = [current;coilcur{i}.'];
    end
    if ismgrid
        numwind = coil_or_mgrid.raw_coil_cur(i);
        coil_or_mgrid.scale_factor(i) = tval;  %Think numwind is already done!
    end

end
if iscoil
    coil_or_mgrid.coilcur = coilcur;
    coil_or_mgrid.current = current;
end

