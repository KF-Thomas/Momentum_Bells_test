function [E_val, corrs] = E(ports,opts_E)
%calculates the quantum correlator E for our Bell test
%input data as structure of the 4 different ports
corr_opts.verbose = opts_E.verbose;
%% options for correlation caculation
if opts_E.plots
    corr_opts.plots=true;
    corr_opts.fig='top halo bb corr';
    corr_opts.direction_labels = {'z','x','y'};
else
    corr_opts.plots=false;
end

corr_opts.type='radial_bb';
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*2e-3;
one_d_range=0.02;
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));
corr_opts.rad_smoothing=nan;

corr_opts.low_mem=true;

corr_opts.sampling_method='complete';
corr_opts.sample_proportion=0.5; %proportion of uncorolated pairs to calculate 
corr_opts.attenuate_counts=1; %artifical qe
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;

corr_opts.timer=false;
corr_opts.print_update = false;

corr_opts.fit = true;


corr_opts.calc_err = opts_E.calc_err;
corr_opts.samp_frac_lims=[0.25,0.5];
corr_opts.num_samp_frac=2;
corr_opts.num_samp_rep=5;

%%
%indexing of ports
%1: top right
%2: top left
%3: bottom right
%4: bottom left

counts12 = [ports.top_right';ports.top_left'];
counts13 = [ports.top_right';ports.bottom_right'];
counts14 = [ports.top_right';ports.bottom_left'];
counts23 = [ports.top_left';ports.bottom_right'];
counts24 = [ports.top_left';ports.bottom_left'];
counts34 = [ports.bottom_right';ports.bottom_left'];

if opts_E.plots
    corr_opts.fig='g14';
end
corrs.g14 = calc_any_g2_type(corr_opts,counts14);
if opts_E.plots
    corr_opts.fig='g23';
end
corrs.g23 = calc_any_g2_type(corr_opts,counts23);
if opts_E.plots
    corr_opts.fig='g12';
end
corrs.g12 = calc_any_g2_type(corr_opts,counts12);
if opts_E.plots
    corr_opts.fig='g34';
end
corrs.g34 = calc_any_g2_type(corr_opts,counts34);
if opts_E.plots
    corr_opts.fig='g13';
end
corrs.g13 = calc_any_g2_type(corr_opts,counts13);
if opts_E.plots
    corr_opts.fig='g24';
end
corrs.g24 = calc_any_g2_type(corr_opts,counts24);

%Caculate the correlator
E_val = (corrs.g14.norm_g2.fitted_g2peak+corrs.g23.norm_g2.fitted_g2peak-corrs.g12.norm_g2.fitted_g2peak-corrs.g34.norm_g2.fitted_g2peak)/...
    (corrs.g14.norm_g2.fitted_g2peak+corrs.g23.norm_g2.fitted_g2peak+corrs.g12.norm_g2.fitted_g2peak+corrs.g34.norm_g2.fitted_g2peak);
end