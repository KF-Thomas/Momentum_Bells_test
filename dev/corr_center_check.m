function out=corr_center_check(counts,fig_title)

out = zeros(3,2);
%% BACK TO BACK (in the same halo)
corr_opts.type='1d_cart_bb';%

corr_opts.verbose = false;
corr_opts.print_update = false;
corr_opts.timer=false;

corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*7e-3;
one_d_range=0.017;%0.02
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
corr_opts.rad_smoothing=nan;
direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

corr_opts.sample_proportion=1.0;%0.65;%1500;
corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')

corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=nan;
corr_opts.sort_norm=0;
corr_opts.one_d_smoothing=0;

if isequal(corr_opts.type,'1d_cart_cl')  || isequal(corr_opts.type,'1d_cart_bb')
    corr_func=@corr_1d_cart;
    corr_density='one_d_corr_density';
    centers='x_centers';
    corr_unc='corr_unc';
    
elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
    corr_func=@corr_radial;
    corr_density='rad_corr_density';
    centers='rad_centers';
    corr_unc='rad_corr_unc';
    direction_label='r';
elseif isequal(corr_opts.type,'3d_cart_cl')  || isequal(corr_opts.type,'3d_cart_bb')
    warning('3d corrs temporarily removed')
end

if isequal(corr_opts.type(end-2:end),'_cl')
    corr_opts.cl_or_bb=false;
elseif isequal(corr_opts.type(end-2:end),'_bb')
    corr_opts.cl_or_bb=true;
end

for ii = 1:3
    
    corr_opts.one_d_dimension = ii;
    direction_label=direction_labels{corr_opts.one_d_dimension};
    
    counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,1);
    normcorr=corr_func(corr_opts,counts_chunked);
    
%     stfig('corr center check')
    stfig(fig_title);
    if ii == 1
        clf
    end
    
    subplot(3,1,ii)
    plot(normcorr.(centers).*1e3,normcorr.(corr_density),'.k-','MarkerSize',10)
    hold on
    fun1d =  @(b,x) b(1).*exp(-((x-b(3)).^2)./(2*b(2).^2))+b(4);
%     fun1d =  @(b,x) b(1)./((x-b(3)).^2+b(2).^2);
    [muHat,sigmaHat] = normfit(normcorr.(centers),0.01,zeros(size(normcorr.(centers))),abs(normcorr.(corr_density)).^2);
    inital_guess=[max(normcorr.(corr_density)),sigmaHat,muHat,0];
%     fit=fitnlm(normcorr.(centers),normcorr.(corr_density),...
%         fun1d,...
%         inital_guess);
    out(ii,1) = trapz(normcorr.(centers),normcorr.(centers).*normcorr.(corr_density))/trapz(normcorr.(centers),normcorr.(corr_density));%fit.Coefficients.Estimate(3);
%     out(ii,2) = fit.Rsquared.Adjusted;
    xx = linspace(-max(normcorr.(centers)),max(normcorr.(centers)),3e3)';
    
%     [ypred,ypredci] = predict(fit,xx,'Simultaneous',true);
%     plot(xx,ypred,'b-', xx,ypredci,'r-');
    
    ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
    xlabel(sprintf('$\\Delta %s$ Seperation $\\times 10^{-3}$',direction_label))
    xlim([-max(normcorr.(centers)),max(normcorr.(centers))].*1e3)
end

end