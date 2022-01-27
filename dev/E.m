function [E_val, corrs] = E(ports,opts_E)
%calculates the quantum correlator E for our Bell test
%input data as structure of the 4 different ports
corr_opts.verbose = logical(opts_E.verbose);
%% general options for correlation caculation
if opts_E.plots
    corr_opts.plots=true;
    corr_opts.fig='place holder';
    corr_opts.direction_labels = {'z','x','y'};
else
    corr_opts.plots=false;
end

%volume widths
dkx = opts_E.delta_kd(2);
dky = opts_E.delta_kd(3);
dkz = opts_E.delta_kd(1);
dkr = (dkx.*dky.*dkz).^(1/3);

% 
corr_density='one_d_corr_density';
if isfield(opts_E,'vol_corr') && opts_E.vol_corr
    corr_opts.type='1d_vol_bb';%'1d_cart_bb';%'2d_cart_bb';%'radial_bb';%
    corr_opts.bin_lims = opts_E.bin_lim;
else
    corr_opts.type='1d_cart_bb';%'radial_bb';
end


% corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]].*0.005;
corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
if isfield(opts_E,'dim')
    corr_opts.one_d_dimension=opts_E.dim;
else
    corr_opts.one_d_dimension=1;%3;%2;
end

if isfield(opts_E,'do_norm')
    corr_opts.g2_norm = opts_E.do_norm;
else
    corr_opts.g2_norm = true;
end
    

if opts_E.norm
    one_d_range=0.16;
    data_type = 'norm';
else
    one_d_range=0.03;%0.015;%0.022;%0.017;%0.06;%
    data_type = 'unnorm';
end


num_pts_cart = round(one_d_range./opts_E.delta_kd(corr_opts.one_d_dimension));
num_pts_rad = round(one_d_range./dkr);
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);%
% one_d_edges_half = [-flip(1:2:num_pts_cart),1:2:num_pts_cart].*opts_E.delta_kd(corr_opts.one_d_dimension);
% corr_opts.one_d_edges = [-flip(1:2:2*num_pts_cart),1:2:2*num_pts_cart].*opts_E.delta_kd(corr_opts.one_d_dimension)./2;%
corr_opts.rad_smoothing=nan;

corr_opts.low_mem=true;

corr_opts.sampling_method='complete';
if isfield(opts_E,'sample_proportion')
    corr_opts.sample_proportion=opts_E.sample_proportion; %proportion of uncorolated pairs to calculate 
else
    corr_opts.sample_proportion=1.0; %proportion of uncorolated pairs to calculate 
end
corr_opts.attenuate_counts=1; %artifical qe
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;

corr_opts.timer=false;
corr_opts.print_update = false;

corr_opts.fit = opts_E.fit;
corr_opts.gaussian_fit = true;

corr_opts.calc_err = opts_E.calc_err;
corr_opts.samp_frac_lims=[0.25,0.5];
corr_opts.num_samp_frac=2;
corr_opts.num_samp_rep=100;%500;

port_pairs = {'g14','g23','g12','g34','g13','g24'};

for this_port = port_pairs
    corr_opts.(this_port{1}) = corr_opts;
    corr_opts.(this_port{1}).fig = this_port{1};
end

%% Port specific options

% corr_opts.g14.redges=sqrt(linspace(0,0.02^2,40));
% corr_opts.g14.rad_smoothing=nan;
% 
% corr_opts.g23.redges=sqrt(linspace(0,one_d_range^2,40));
% corr_opts.g23.rad_smoothing=3e-5;

%%
%indexing of ports
%1: top right
%2: top left
%3: bottom right
%4: bottom left

counts12 = [ports.top_right.(data_type)';ports.top_left.(data_type)'];
counts13 = [ports.top_right.(data_type)';ports.bottom_right.(data_type)'];
counts14 = [ports.top_right.(data_type)';ports.bottom_left.(data_type)'];
counts23 = [ports.top_left.(data_type)';ports.bottom_right.(data_type)'];
counts24 = [ports.top_left.(data_type)';ports.bottom_left.(data_type)'];
counts34 = [ports.bottom_right.(data_type)';ports.bottom_left.(data_type)'];

corrs.g14 = calc_any_g2_type(corr_opts.g14,counts14); %btw halo bb

corrs.g23 = calc_any_g2_type(corr_opts.g23,counts23); %btw halo bb

corrs.g12 = calc_any_g2_type(corr_opts.g12,counts12); %top halo bb

corrs.g34 = calc_any_g2_type(corr_opts.g34,counts34); %btm halo bb

% corrs.g13 = calc_any_g2_type(corr_opts.g13,counts13);

% corrs.g24 = calc_any_g2_type(corr_opts.g24,counts24);

%Caculate the correlator
if corr_opts.g2_norm
mid_pt=ceil(length(corrs.g12.norm_g2.g2_amp)/2);
g12 = corrs.g12.norm_g2.g2_amp(mid_pt);
g14 = corrs.g14.norm_g2.g2_amp(mid_pt);
g23 = corrs.g23.norm_g2.g2_amp(mid_pt);
g34 = corrs.g34.norm_g2.g2_amp(mid_pt);
else
    mid_pt=ceil(length(corrs.g12.in_shot_corr.(corr_density))/2);
g12 = corrs.g12.in_shot_corr.(corr_density)(mid_pt);
g14 = corrs.g14.in_shot_corr.(corr_density)(mid_pt);
g23 = corrs.g23.in_shot_corr.(corr_density)(mid_pt);
g34 = corrs.g34.in_shot_corr.(corr_density)(mid_pt);
end

E_val = (g14+g23-g12-g34)/...
    (g14+g23+g12+g34);
end