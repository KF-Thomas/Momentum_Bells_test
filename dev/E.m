function [E_val, corrs] = E(ports,opts_E)
%calculates the quantum correlator E for our Bell test
%input data as structure of the 4 different ports

%% options for correlation caculation
if opts_E.plots
corr_opts.plots=true;
corr_opts.fig='top halo bb corr';
corr_opts.direction_labels = {'z','x','y'};
end

corr_opts.type='radial_bb';
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*2e-3;
one_d_range=0.02;
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));
corr_opts.rad_smoothing=nan;

corr_opts.low_mem=true;

corr_opts.norm_samp_factor=1500;
corr_opts.attenuate_counts=1;
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
counts14 = [ports.top_right';ports.bottom_left'];
counts23 = [ports.top_left';ports.bottom_right'];
counts34 = [ports.bottom_right';ports.bottom_left'];


corrs.g14 = calc_any_g2_type(corr_opts,counts14);
corrs.g23 = calc_any_g2_type(corr_opts,counts23);
corrs.g12 = calc_any_g2_type(corr_opts,counts12);
corrs.g34 = calc_any_g2_type(corr_opts,counts34);


b = fit.Coefficients.Estimate;
corrs.top_halo.corr_bb.fit
E_val = ;
end