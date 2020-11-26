function out=center_check(counts)
%% BACK TO BACK (in the same halo)
corr_opts.type='1d_cart_bb';%

corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*3e-3;
one_d_range=0.017;%0.02
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
corr_opts.rad_smoothing=nan;
direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

corr_opts.sample_proportion=1;%0.65;%1500;
corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')

corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=nan;
corr_opts.sort_norm=0;

for ii = 1:3

corr_opts.one_d_dimension = ii;

if isequal(corr_opts.type,'1d_cart_cl')  || isequal(corr_opts.type,'1d_cart_bb')
    corr_func=@corr_1d_cart;
    corr_density='one_d_corr_density';
    centers='x_centers';
    corr_unc='corr_unc';
    direction_label=direction_labels{corr_opts.one_d_dimension};
    bin_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2;
elseif isequal(corr_opts.type,'radial_cl')  || isequal(corr_opts.type,'radial_bb')
    corr_func=@corr_radial;
    corr_density='rad_corr_density';
    centers='rad_centers';
    corr_unc='rad_corr_unc';
    direction_label='r';
    bin_centers=sqrt((corr_opts.redges(2:end).^3-corr_opts.redges(1:end-1).^3)./(3.*(corr_opts.redges(2:end)-corr_opts.redges(1:end-1))));
elseif isequal(corr_opts.type,'3d_cart_cl')  || isequal(corr_opts.type,'3d_cart_bb')
    warning('3d corrs temporarily removed')
end


counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,norm_sort_dir);
normcorr=corr_func(corr_opts,counts_chunked);
subplot(3,1,ii)
plot(normcorr.(centers),normcorr.(corr_density),'.k-','MarkerSize',10)

ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
    xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
end

end