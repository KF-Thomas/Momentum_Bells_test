function [E_val, corrs] = E(ports,opts_E)
%calculates the quantum correlator E for our Bell test
%input data as structure of the 4 different ports
corr_opts.verbose = opts_E.verbose;
%% general options for correlation caculation
if opts_E.plots
    corr_opts.plots=true;
    corr_opts.fig='top halo bb corr';
    corr_opts.direction_labels = {'z','x','y'};
else
    corr_opts.plots=false;
end

corr_opts.type='radial_bb';
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*2e-3;
one_d_range=0.017;
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));
corr_opts.rad_smoothing=nan;

corr_opts.low_mem=true;

corr_opts.sampling_method='complete';
corr_opts.sample_proportion=1; %proportion of uncorolated pairs to calculate 
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

port_pairs = {'g14','g23','g12','g34','g13','g24'};

for this_port = port_pairs
    corr_opts.(this_port{1}) = corr_opts;
    corr_opts.(this_port{1}).fig = this_port{1};
end

%% Port specific options

corr_opts.g14.redges=sqrt(linspace((5e-3)^2,0.02^2,40));
corr_opts.g14.rad_smoothing=nan;

corr_opts.g23.redges=sqrt(linspace((5e-3)^2,one_d_range^2,40));
corr_opts.g23.rad_smoothing=3e-5;

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

corrs.g14 = calc_any_g2_type(corr_opts.g14,counts14);

corrs.g23 = calc_any_g2_type(corr_opts.g23,counts23);

corrs.g12 = calc_any_g2_type(corr_opts,counts12);

corrs.g34 = calc_any_g2_type(corr_opts,counts34);

corrs.g13 = calc_any_g2_type(corr_opts.g13,counts13);

corrs.g24 = calc_any_g2_type(corr_opts.g24,counts24);

%Caculate the correlator
E_val = (corrs.g14.norm_g2.fitted_g2peak+corrs.g23.norm_g2.fitted_g2peak-corrs.g12.norm_g2.fitted_g2peak-corrs.g34.norm_g2.fitted_g2peak)/...
    (corrs.g14.norm_g2.fitted_g2peak+corrs.g23.norm_g2.fitted_g2peak+corrs.g12.norm_g2.fitted_g2peak+corrs.g34.norm_g2.fitted_g2peak);
end