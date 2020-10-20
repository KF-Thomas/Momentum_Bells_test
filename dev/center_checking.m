% cartesian g2 normalisation

corr_opts.one_d_dimension = 2;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*3e-3;
one_d_range=0.017;%0.02
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
corr_opts.rad_smoothing=nan;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

    corr_func=@corr_1d_cart;
    corr_density='one_d_corr_density';
    centers='x_centers';
    corr_unc='corr_unc';
    direction_label=corr_opts.direction_labels{corr_opts.one_d_dimension};
%     bin_centers=(corr_opts.redges(2:end)+corr_opts.redges(1:end-1))/2;

counts = top_halo.counts_vel';
counts = bottom_halo.counts_vel';

counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,norm_sort_dir);
normcorr=corr_func(corr_opts,counts_chunked);