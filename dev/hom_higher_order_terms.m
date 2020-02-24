

% HOM x halo
% the probability of detecting k particles with lambda mean occupation is
p_good=@(mode_occ) poisson_dist(1,mode_occ).*poisson_dist(1,mode_occ);


mode_samp=logspace(-2.5,0,1e2);

p_ratio=@(mode_occ)   p_good(mode_occ)./p_bad(mode_occ);

p_ratio_samp=p_ratio(mode_samp);

figure(1)
loglog(mode_samp,p_ratio_samp)
xlabel('mode occ')
ylabel('prob good against bad')

% the better way to do it is to find the chance of a coincidence across all output states of the interferometer
% P_CO = 1- P(|0,N>)- P(|N,0>)
% with some atom number in port a and b

%%
figure(2)
nmax=30; % increase to consider larger atom totals
clf
set(gcf,'color','w')
loglog(mode_samp,mode_samp,'g-')
hold on
coherent_co=coincidence_coherent(mode_samp,nmax);
incoherent_co=coincidence_incoherent(mode_samp,nmax);
loglog(mode_samp,coherent_co,'k-')
loglog(mode_samp,incoherent_co,'b-')
hold off
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('mode occ')
ylabel('coincidence probability')
legend('y=x','coherent','incoherent','Location','northwest')

%%
coincidence_coherent([0.1])
coincidence_coherent([0.4])
coincidence_coherent([0.1,0.2,0.3,0.4])
coincidence_weighted([0.1,0.2,0.3,0.4])


%% we can also have a look at how quantum supresses coincidences (cf classical non interacting) for higher 
% order HOM
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

% plots
coinc.ntot=coinc.na+coinc.nb;
zvalue=abs(coinc.na-coinc.nb)./coinc.ntot;
zvalue_scaled=zvalue-min(zvalue);
zvalue_scaled=zvalue_scaled/max(zvalue_scaled);
cmap=parula(1e3);
zvalue_scaled=1+zvalue_scaled*(size(cmap,1)-2);
figure(3)
clf
scatter(coinc.ntot,coinc.quantum,30, cmap(ceil(zvalue_scaled),:), 'filled','o','MarkerEdgeColor',[0,0,0])
hold on
plot(coinc.ntot,coinc.class,'-')
plot(coinc.quant_avg.ntot,coinc.quant_avg.coinc,'-')
legend('quantum','classical particles','average for quantum with n total')
hold off
xlabel('total number in')
ylabel('coincidence probability')
%set(gca,'XScale','log')
%set(gca,'YScale','log')
hcb=colorbar
hcb.Label.String = 'input number factional disparity 0=equal, 1=all one port';



%%
coincidence_coherent(0.1,1)
%%
% some prior art but its not exactly what im after
%https://iopscience.iop.org/article/10.1088/1367-2630/ab1bbf

function state=out_state_beam_splitter(m,n,distinguishible)
    if nargin<3
        distinguishible=0;
    end
    % given an inital state (a_dag)^n (b_dag)^m |0,0>_{a,b}
    % where a_dag and b_dag are the creation operators in mode a and b
    % to find the output state we do the expansion (c_dag+d_dag)^n (c_dag-d_dag)^m |0,0>_{c,d}
    % where c_dag and d_dag are the creation operators in mode c and d 
    % we will buld the output state in a matrix representation
    % where the first index is the power of c +1 and the second is the power of d +1, 
    % this is bc of index starting at 1 and we need to rep. 0
    % the value in the matrix is the prefactor
    % the double loop is a way of doing the distributive multipication
    coef_amp_ab=zeros(m+n+1);
    for ii=0:m
        c_pow_a=m-ii;
        d_pow_a=ii;
        % find the amplitude using the binomial theorem
        coef_amp_a=nchoosek(m,ii);
        for jj=0:n
            c_pow_b=n-jj;
            d_pow_b=jj;
            %https://math.stackexchange.com/questions/1861168/is-there-a-formula-for-the-binomial-expansion-of-a-bn
            if distinguishible
                coef_amp_b=nchoosek(n,jj);
            else
                coef_amp_b=((-1)^(jj))*nchoosek(n,jj);
            end
            %fprintf('A side  %.0f c^%u d^%u \n',coef_amp_a,c_pow_a,d_pow_a)
            %fprintf('B side  %.0f c^%u d^%u \n',coef_amp_b,c_pow_b,d_pow_b)
            %fprintf('element %.0f c^%u d^%u \n',coef_amp_b*coef_amp_a,c_pow_b+c_pow_a,d_pow_b+d_pow_a)
            %pause
            %fprintf('----\n')
            % the position in the matrix is the sum of the powers
            % and the value is the product of the coefs added on to what was already in the matrix
            coef_amp_ab(c_pow_a+c_pow_b+1,d_pow_a+d_pow_b+1)=coef_amp_ab(c_pow_a+c_pow_b+1,d_pow_a+d_pow_b+1)+...
                                                                coef_amp_a*coef_amp_b;
        end
    end
    % the sum of amplituded is then
    coef_amp=abs(coef_amp_ab);
    total_amp=sum(coef_amp(:));
    % the coincidence probability is then the total amplitude minus  P(|0,N>) and P(|N,0>) 
    % (no atoms in the other port)
    state=coef_amp_ab;
   
end



%%

% easier to calculate coincidence than g2
function co_out=coincidence_coherent(mode_occ,n_max)
if nargin<2
    n_max=50;
end
% lets sum over all ocupations
% the factorial accounts for the colinear correlation increasing the chance of getting n particles
% (over random chance) due to thermal bunching
% i think this approach misses the interference of output states
state_out_tot=zeros(n_max+1,n_max+1,numel(mode_occ));
for n_tot=1:n_max
   for num_b=0:n_tot
       num_a=n_tot-num_b;
       % the statistical chance of getting this number combination
       stat_prob_nums=poisson_dist(num_a,mode_occ).*poisson_dist(num_b,mode_occ);
       % the quantum enhancement of this combination
       qunat_prob_nums=factorial(num_a).*factorial(num_b).*stat_prob_nums;
       % the coincidence probability from HOM given this number combination in
       fprintf('num a %u, num b %u \n',num_a,num_b);
       component_state=out_state_beam_splitter(num_a,num_b)
       for kk=1:numel(mode_occ) %step over the query mode occupancies
            % padd the beamsplitter output state weight it by qunat_prob_nums and add it to the output state
            state_out_tot(:,:,kk)=state_out_tot(:,:,kk)+padarray(component_state*qunat_prob_nums(kk),[1,1]*(n_max+1-(num_a+num_b+1)),0,'post');
       end
   end
end
state_out_tot
out_state_amp=abs(state_out_tot);
co_out=( sum(out_state_amp,[1,2])- sum(out_state_amp(:,1,:),1)+sum(out_state_amp(1,:,:),2) )./sum(out_state_amp,[1,2]);
co_out=squeeze(co_out);
end

% easier to calculate coincidence than g2
function co_out=coincidence_incoherent(mode_occ,n_max)
if nargin<2
    n_max=50;
end
% lets sum over all ocupations
% the factorial accounts for the colinear correlation increasing the chance of getting n particles
% (over random chance) due to thermal bunching
% i think this approach misses the interference of output states
out_state=zeros(n_max+1,n_max+1,numel(mode_occ));
for ii=1:n_max
   for jj=0:ii
       num_a=ii-jj;
       num_b=jj;
       % the statistical chance of getting this number combination
       stat_prob_nums=poisson_dist(num_a,mode_occ).*poisson_dist(num_b,mode_occ);
       % the quantum enhancement of this combination
       qunat_prob_nums=factorial(num_a).*factorial(num_b).*stat_prob_nums;
       % the coincidence probability from HOM given this number combination in
       [~,state]=coincidence_prob_incoherent(num_a,num_b);
       for kk=1:numel(mode_occ) %step over the query mode occupancies
            % padd the beamsplitter output state weight it by qunat_prob_nums and add it to the output state
            out_state(:,:,kk)=out_state(:,:,kk)+padarray(state*qunat_prob_nums(kk),[1,1]*(n_max+1-(num_a+num_b+1)),0,'post');
       end
   end
end
out_state_amp=abs(out_state);
co_out=( sum(out_state_amp,[1,2])- sum(out_state_amp(:,1,:),1)+sum(out_state_amp(1,:,:),2) )./sum(out_state_amp,[1,2]);
co_out=squeeze(co_out);
end

%%

% easier to calculate coincidence than g2
function co_out=coincidence_weighted(mode_occ,n_max)
if nargin<2
    n_max=50;
end
n_max=50;
% lets sum over all ocupations
% the factorial accounts for the colinear correlation increasing the chance of getting n particles
% (over random chance) due to thermal bunching
% i think this approach misses the interference of output states
co_out=0;
for ii=1:n_max
   for jj=0:ii
       num_a=ii-jj;
       num_b=jj;
       % the statistical chance of getting this number combination
       stat_prob_nums=poisson_dist(num_a,mode_occ).*poisson_dist(num_b,mode_occ);
       % the quantum enhancement of this combination
       qunat_prob_nums=factorial(num_a).*factorial(num_b).*stat_prob_nums;
       % the coincidence probability from HOM given this number combination in
       hom_coincidence=coincidence_prob(num_a,num_b);
       co_out=co_out+qunat_prob_nums*hom_coincidence;
   end
end


end



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






