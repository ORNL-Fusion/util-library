function coil = set_w7x_current(coil,taper)
           %1  2  3  4  5  6  7   8   9       10       11       12       13       14
% taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2,Itrim_A1,Itrim_A2,Itrim_A3,Itrim_A4,Itrim_B1];
% taper should be the current per winding, as it is multiplied by the number of windings in each coil belowuw582
taper_tmp = taper;
taper = zeros(1,14);
taper(1:length(taper_tmp)) = taper_tmp;

nc = length(coil.coilpos);
current = [];
icount = 0;
for i = 1:nc
    cname = lower(char(coil.cname{i}));
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
    coilcur{i} = tval*ones(size(coil.coilcur{i}))*coil.numwind(i);   
    coilcur{i}(end) = 0;
    current = [current;coilcur{i}.'];
    icount = icount + length(coil.coilcur{i});
end

coil.coilcur = coilcur;
coil.current = current;

