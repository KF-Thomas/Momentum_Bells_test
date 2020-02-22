

% HOM x halo
% the probability of detecting k particles with lambda mean occupation is


p_good=@(mode_occ) poisson_dist(1,mode_occ).*poisson_dist(1,mode_occ);


mode_samp=logspace(-4,0,1e3);

p_ratio=@(mode_occ)   p_good(mode_occ)./p_bad(mode_occ);

p_ratio_samp=p_ratio(mode_samp);

figure(1)
loglog(mode_samp,p_ratio_samp)
xlabel('mode occ')
ylabel('prob good against bad')



%
figure(2)
loglog(mode_samp,coincidence_weighted(mode_samp))
xlabel('mode occ')
ylabel('coincidence')


%%
coincidence_prob(1,1)
%% look at how quantum supresses coincidences
coinc={};
coinc.quantum=[];
coinc.class=[];
coinc.na=[];
coinc.nb=[];
ii=1;
for ntot=1:20
    for nb=0:ntot
        na=ntot-nb;
        coinc.quantum(ii)=coincidence_prob(na,nb);
        coinc.class(ii)=1-(1/2)^ntot;
        coinc.na(ii)=na;
        coinc.nb(ii)=nb;
        ii=ii+1;
	end
end

% find mean conincidence for quantum case with ntot
coinc.ntot=coinc.na+coinc.nb;
nvals=unique(coinc.ntot);
coinc.quant_avg=[];
coinc.quant_avg.ntot=nvals;
coinc.quant_avg.coinc=nvals*nan;
for ii=1:numel(nvals)
    coinc.quant_avg.coinc(ii)=mean(coinc.quantum(coinc.ntot==(nvals(ii))));
end

%%
coinc.ntot=coinc.na+coinc.nb;
zvalue=abs(coinc.na-coinc.nb)./coinc.ntot;
zvalue_scaled=zvalue-min(zvalue);
zvalue_scaled=zvalue_scaled/max(zvalue_scaled);
cmap=parula(1e3);
zvalue_scaled=1+zvalue_scaled*(size(cmap,1)-2);
figure(1)
clf
scatter(coinc.ntot,coinc.quantum,30, cmap(ceil(zvalue_scaled),:), 'filled','o','MarkerEdgeColor',[0,0,0])
hold on
plot(coinc.ntot,coinc.class,'-')
plot(coinc.quant_avg.ntot,coinc.quant_avg.coinc,'-')
legend('quantum','classical particles','average for quantum with n total')
hold off
%set(gca,'XScale','log')
%set(gca,'YScale','log')
hcb=colorbar
hcb.Label.String = 'input number factional disparity 0=equal, 1=all one port';



%%
coincidence_prob(3,3)
%%

%https://iopscience.iop.org/article/10.1088/1367-2630/ab1bbf
% i cant get that exactly to work
function prob=coincidence_prob(m,n)
% single photon coincidence
if m==0 && n==0
    prob=0;
else
   coef_amp_ab=zeros(m+n+1);

    % this basicaly just does the expansion (c+d)^n (c-d)^m
    % and represents the expanded state in a matrix representation where the first index is the power of c +1
    % and the second is the power of d +1, the value in the matris is the prefactor
    % the double loop is a way of doing the distributive multipication
    for ii=0:m
        c_pow_a=m-ii;
        d_pow_a=ii;
        % find the amplitude using the binomial theorem
        coef_amp_a=nchoosek(m,ii);
        for jj=0:n
            c_pow_b=n-jj;
            d_pow_b=jj;
            %https://math.stackexchange.com/questions/1861168/is-there-a-formula-for-the-binomial-expansion-of-a-bn
            coef_amp_b=((-1)^(jj))*nchoosek(n,jj);
            %fprintf('A side  %.0f c^%u d^%u \n',coef_amp_a,c_pow_a,d_pow_a)
            %fprintf('B side  %.0f c^%u d^%u \n',coef_amp_b,c_pow_b,d_pow_b)
            %fprintf('element %.0f c^%u d^%u \n',coef_amp_b*coef_amp_a,c_pow_b+c_pow_a,d_pow_b+d_pow_a)
            %pause
            %fprintf('----\n')
            coef_amp_ab(c_pow_a+c_pow_b+1,d_pow_a+d_pow_b+1)=coef_amp_ab(c_pow_a+c_pow_b+1,d_pow_a+d_pow_b+1)+...
                                                                coef_amp_a*coef_amp_b;
        end
    end
    coef_amp=abs(coef_amp_ab);
    total_amp=sum(coef_amp(:));
    prob=(total_amp-sum(coef_amp(:,1))-sum(coef_amp(1,:)))/total_amp;
end
   
end


%%

function p_out=poisson_dist(k,lambda)
    p_out=(lambda.^k).*exp(-lambda)./ factorial(k);
end




function p_out=p_bad(mode_occ)
n_max=100;
% find the probability of having 2 in one mode
% this will be double posonian bc of g^2=2
% now bc this wont contribute to the x halo g2 i dont think it should be included
% p_out=2*poisson_dist(2,mode_occ).*poisson_dist(0,mode_occ)+...
%         2*poisson_dist(0,mode_occ).*poisson_dist(2,mode_occ);
p_out=0;

for ii=3:n_max
   for jj=0:ii
       p_out=p_out+factorial(ii-jj)*poisson_dist(ii-jj,mode_occ).*factorial(jj).*poisson_dist(jj,mode_occ);
   end
end


end

% easier to calculate coincidence than g2
function co_out=coincidence_weighted(mode_occ)
n_max=50;
% lets sum over all ocupations
% the factorial accounts for the colinear correlation increasing the chance of getting n particles
% due to thermal bunching
co_out=0;
for ii=1:n_max
   for jj=0:ii
       num_a=ii-jj;
       num_b=jj;
       stat_prob_nums=poisson_dist(num_a,mode_occ).*poisson_dist(num_b,mode_occ);
       qunat_prob_nums=factorial(num_a).*factorial(num_b).*stat_prob_nums;
       hom_coincidence=coincidence_prob(num_a,num_b);
       co_out=co_out+qunat_prob_nums*hom_coincidence;
   end
end


end




