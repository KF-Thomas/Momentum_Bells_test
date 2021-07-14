function corrs = global_corrs(top_halo,bottom_halo,global_opts) 
%% general settings
corr_opts.verbose = false;
corr_opts.print_update = false;
corr_opts.timer=false;

corr_opts.plots = global_opts.plots;
corr_opts.fit = global_opts.fit;
corr_opts.calc_err = global_opts.calc_err;

global_sample_portion = 1.0;%0.05;%0.08;%0.5;%1.0;%0.001;%

% variables for calculating the error
corr_opts.samp_frac_lims=[0.65,0.9];
corr_opts.num_samp_frac=5;
corr_opts.num_samp_rep=5;

corr_opts.attenuate_counts=1;

%% BACK TO BACK (in the same halo)
corr_opts.type='1d_cart_bb';%'radial_bb';%
corr_opts.one_d_dimension = 2;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*40e-3;
one_d_range=0.01;%0.017;%0.09;%0.075;%0.02%0.03
% one_d_range=0.16;
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));%100 or 80 or 85 or 95
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';
corr_opts.rad_smoothing=nan;%0.0003;%
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

corr_opts.norm_samp_factor=1500;%1500;
corr_opts.sample_proportion=global_sample_portion;%0.65;%1500;
corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=nan;
corr_opts.sort_norm=0;

corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit

%% TOP HALO BACK TO BACK

corr_opts.fig='top halo bb corr';
corrs.top_halo.corr_bb=calc_any_g2_type(corr_opts,top_halo.counts_vel');
% corrs.top_halo.corr_bb=calc_any_g2_type(corr_opts,top_halo.counts_vel_norm');

%% BOTTOM HALO BACK TO BACK

corr_opts.fig='bottom halo bb corr';
corrs.bottom_halo.corr_bb=calc_any_g2_type(corr_opts,bottom_halo.counts_vel');
% corrs.bottom_halo.corr_bb=calc_any_g2_type(corr_opts,bottom_halo.counts_vel_norm');

%% CO-LINEAR (in the same halo)

corr_opts.type='radial_cl';%'1d_cart_cl';%'3d_cart_cl';%%%
% corr_opts.one_d_dimension=3;
% corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*2e-3;
one_d_range=0.02;
% corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,74)';
corr_opts.redges=sqrt(linspace(0,one_d_range^2,30));
corr_opts.rad_smoothing=nan;
% corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=nan;

corr_opts.sampling_method='complete';%
corr_opts.norm_samp_factor=1500;
corr_opts.sample_proportion=global_sample_portion;%1500;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;

% corr_opts.one_d_smoothing=nan;
% corr_opts.one_d_smoothing=0.002;

%% TOP HALO CO-LINEAR

corr_opts.fig='top halo cl corr';
% corrs.top_halo.corr_cl=calc_any_g2_type(corr_opts,top_halo.counts_vel');

%% BOTTOM HALO CO-LINEAR

corr_opts.fig='bottom halo cl corr';
% corrs.bottom_halo.corr_cl=calc_any_g2_type(corr_opts,bottom_halo.counts_vel');

%% BETWEEN HALO CO-LINEAR
corr_opts.fig='between halo cl corr';
corr_opts.type='radial_cl';%'3d_cart_cl';%'1d_cart_cl';%
corr_opts.one_d_dimension=2;
% corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*0.92;
% one_d_range=0.16;%0.04;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2;
one_d_range=0.02;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,30)';
corr_opts.redges=sqrt(linspace(1e-6^2,one_d_range^2,30));
corr_opts.rad_smoothing=nan;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;
corr_opts.sampling_method='complete';%
corr_opts.sample_proportion=global_sample_portion;%1500;
corr_opts.norm_samp_factor=1;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;
corr_opts.progress_updates=5;

% both_halo_counts = [top_halo.counts_vel_norm';bottom_halo.counts_vel_norm'];
both_halo_counts = [top_halo.counts_vel';bottom_halo.counts_vel'];

corr_opts.one_d_smoothing=nan;
% corr_opts.one_d_smoothing=0.01;

% corrs.between_halos.corr_cl=calc_any_g2_type(corr_opts,both_halo_counts);

%% BETWEEN HALO BACK TO BACK
corr_opts.fig='between halo bb corr';
corr_opts.type='radial_bb';%'3d_cart_bb';%'1d_cart_bb';%
corr_opts.one_d_dimension=3;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*0.05;
% one_d_range=0.16;
% corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2;
one_d_range=0.017;%0.09;%
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,100)';
corr_opts.redges=sqrt(linspace(0,one_d_range^2,75));%(5e-3)^2
corr_opts.rad_smoothing=nan;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;
corr_opts.sampling_method='complete';%'basic';%'complete';%
corr_opts.norm_samp_factor=1500;%1500;
corr_opts.sample_proportion=global_sample_portion;%1500;
corr_opts.norm_samp_factor=1;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;
corr_opts.progress_updates=5;

corr_opts.gaussian_fit = false; %ensure it always uses a gaussian fit

% both_halo_counts = [top_halo.counts_vel_norm';bottom_halo.counts_vel_norm'];
both_halo_counts = [top_halo.counts_vel';bottom_halo.counts_vel'];

corr_opts.one_d_smoothing=nan;
% corr_opts.one_d_smoothing=0.0008;

corrs.between_halos.corr_bb=calc_any_g2_type(corr_opts,both_halo_counts);
end