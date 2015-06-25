clearvars;

verbose = 1;

gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';
g=readg_g3d(gfile_name);
psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);

q = g.qpsi;
pn = g.pn;


nmax = 2;
mmax = 20;


qmin = min(q);
qmax = max(q);
fprintf('\n\n  q profile [min,max] = [%8.2f,%8.2f]\n',qmin,qmax)

figure; hold on; box on;
plot(pn,q,'linewidth',2)

ylim = get(gca,'ylim');
ylab = (ylim(2)-ylim(1))*0.9 + ylim(1);

pn_res_mat_mn = zeros(mmax,nmax);

for m = 1:mmax
    for n = 1:nmax
        qres = m/n;
        if verbose
            fprintf('Testing m = %2i, n=%2i, qres=%12.4f --> ',m,n,qres)
        end
        if (qres < qmax) && (qres > qmin)
            pn_res = interp1(q,pn,qres);
            if ~any(any(pn_res_mat_mn == pn_res))
                text(pn_res,ylab,[num2str(m),'/',num2str(n)])
            end
            pn_res_mat_mn(m,n) = pn_res;
            plot([1,1]*pn_res,ylim,'k')                            
            
            if verbose
                fprintf(' res at pn = %12.4f\n',pn_res)            
            end
        else
            if verbose
                fprintf(' NO\n')
            end
        end
        if (qres < qmin) 
            break;
        end
    end
end

title('m/n resonances','fontsize',12)
xlabel('\psi_N','fontsize',12)            
ylabel('q','fontsize',12)            
set(gca,'fontsize',12)