clearvars;

geo_path = 'C:/Work/Stellarator/All_W7X_WORK/geometry_files/';
mypath = 'C:\Work_archive\RUN_ARCHIVE\DIV3D_RUNS\W7X\OP1.2a\22kA_9.155e-6\';

% thispart_1 = [1,3,9,21];  %"Everything -- remove occluding baffles, no scraper"
% thispart_5 = [5,7,11,13,19];  %

thispart_1 = [1,3,21];  %Just basics
thispart_5 = [19];  %

nex_tor = ones(1,50)*3;
% nex_pol = zeros(1,50);
nex_pol = ones(1,50)*5;
nex_pol([19,21]) = 10*5;

stellsym = 1;
nfp = 5;


FP_WANT1 = 3;
MOVE_1 = 1;
FP_WANT2 = 4;
MOVE_2 = 0;

plotit = 0;
newfig = 0;
[xpart,ypart,zpart,~,~,ntor,npol]=load_allparts(mypath,stellsym,[],plotit,newfig,[],[],nfp);

% Simple strip of last poloidal on each part  - skip one side of 1
if 1
    for k = 1:size(xpart,1)
        if k == 1
            xpart(k,1:ntor(k),1:npol(k)-1) = xpart(k,1:ntor(k),2:npol(k));
            ypart(k,1:ntor(k),1:npol(k)-1) = ypart(k,1:ntor(k),2:npol(k));
            zpart(k,1:ntor(k),1:npol(k)-1) = zpart(k,1:ntor(k),2:npol(k));
            npol(k) = npol(k) - 1;
        else
            xpart(k,1:ntor(k),1:npol(k)-2) = xpart(k,1:ntor(k),2:npol(k)-1);
            ypart(k,1:ntor(k),1:npol(k)-2) = ypart(k,1:ntor(k),2:npol(k)-1);
            zpart(k,1:ntor(k),1:npol(k)-2) = zpart(k,1:ntor(k),2:npol(k)-1);
            npol(k) = npol(k) - 2;
        end
    end
end

% Now make a small offset
if 1
DEBUG = 0;
dX = 0.00001;
if DEBUG
figure; hold on;
end

for k = 1:size(xpart,1)
    
    if DEBUG
        plot3(squeeze(xpart(k,1:ntor(k),1:npol(k)))  ,squeeze(ypart(k,1:ntor(k),1:npol(k)))  ,squeeze(zpart(k,1:ntor(k),1:npol(k)))  ,'k-')
        plot3(squeeze(xpart(k,1:ntor(k),1:npol(k))).',squeeze(ypart(k,1:ntor(k),1:npol(k))).',squeeze(zpart(k,1:ntor(k),1:npol(k))).','k-')
    end
    for i = 1:ntor(k)
        for j = 1:npol(k)
            p1 = [xpart(k,i,j),ypart(k,i,j),zpart(k,i,j)];
            s1 = 1;
            s2 = 1;
            if i < ntor(k)
                p2 = [xpart(k,i+1,j),ypart(k,i+1,j),zpart(k,i+1,j)];
            else
                p2 = [xpart(k,i-1,j),ypart(k,i-1,j),zpart(k,i-1,j)];
                s1 = -1;
            end
            if j < npol(k)
                p3 = [xpart(k,i,j+1),ypart(k,i,j+1),zpart(k,i,j+1)];
            else
                p3 = [xpart(k,i,j-1),ypart(k,i,j-1),zpart(k,i,j-1)];
                s2 = -1;
            end
            
            v12 = p2 - p1;
            v13 = p3 - p1;
            vn = cross(v12,v13);
            p4 = p1 + s1*s2*dX * vn./norm(vn);
            if DEBUG
%                 plot3(p1(1),p1(2),p1(3),'o')
%                 plot3(p2(1),p2(2),p2(3),'x')
%                 plot3(p3(1),p3(2),p3(3),'<')
                plot3(p4(1),p4(2),p4(3),'d')
            end
            xpart(k,i,j) = p4(1);
            ypart(k,i,j) = p4(2);
            zpart(k,i,j) = p4(3);
            
            
        end
    end
    
end

end

[xpart,ypart,zpart,npol,ntor,xel,yel,zel,area,atri] = refine_part(xpart,ypart,zpart,npol,ntor,nex_pol,nex_tor);
rpart = sqrt(xpart.^2 + ypart.^2);
phipart = atan2(ypart,xpart);
phipart(19,ntor(19),1:npol(19))  = 1e-6;


fp_want = FP_WANT1;
move_to_2nd_hfp = MOVE_1;
[Rs1,Zs1,Ps1] = symmetrize_and_move_to_hfp(rpart,zpart,phipart,nfp,stellsym,fp_want,move_to_2nd_hfp);
Xs1 = Rs1.*cos(Ps1);
Ys1 = Rs1.*sin(Ps1);
figure; hold on; box on;
for k = thispart_1
    plot3(squeeze(Xs1(k,1:ntor(k),1:npol(k)))  ,squeeze(Ys1(k,1:ntor(k),1:npol(k)))  ,squeeze(Zs1(k,1:ntor(k),1:npol(k)))  ,'k-')
    plot3(squeeze(Xs1(k,1:ntor(k),1:npol(k))).',squeeze(Ys1(k,1:ntor(k),1:npol(k))).',squeeze(Zs1(k,1:ntor(k),1:npol(k))).','k-')
end

% fp_want = 5;
% move_to_2nd_hfp = 1;
fp_want = FP_WANT2;
move_to_2nd_hfp = MOVE_2;
[Rs5,Zs5,Ps5] = symmetrize_and_move_to_hfp(rpart,zpart,phipart,nfp,stellsym,fp_want,move_to_2nd_hfp);
Xs5 = Rs5.*cos(Ps5);
Ys5 = Rs5.*sin(Ps5);
for k = thispart_5    
    plot3(squeeze(Xs5(k,1:ntor(k),1:npol(k)))  ,squeeze(Ys5(k,1:ntor(k),1:npol(k)))  ,squeeze(Zs5(k,1:ntor(k),1:npol(k)))  ,'k-')
    plot3(squeeze(Xs5(k,1:ntor(k),1:npol(k))).',squeeze(Ys5(k,1:ntor(k),1:npol(k))).',squeeze(Zs5(k,1:ntor(k),1:npol(k))).','k-')
end

npts = 0;
for k = [thispart_1,thispart_5]
    npts = npts + ntor(k)*npol(k);
end

% Write grid file
outfile = 'C:\Work\FLARE\w7x\grid_vf.dat';
fid = fopen(outfile,'w');
fprintf(fid,'# grid_id = 110\n');
fprintf(fid,'# total resolution:  n_xyz   =  %i\n',npts);
for k = thispart_1
    for i = 1:ntor(k)
        for j = 1:npol(k)
            fprintf(fid,'%18.12f %18.12f %18.12f\n',100*Xs1(k,i,j),100*Ys1(k,i,j),100*Zs1(k,i,j));
        end
    end
end
for k = thispart_5
    for i = 1:ntor(k)
        for j = 1:npol(k)
            fprintf(fid,'%18.12f %18.12f %18.12f\n',100*Xs5(k,i,j),100*Ys5(k,i,j),100*Zs5(k,i,j));
        end
    end
end
fclose(fid);


