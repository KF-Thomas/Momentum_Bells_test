% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '20191115_halos_attempt_3';
% data_folder = '20191104_halos_attempt_1';
opts.import.dir = fullfile(opts.data_root, data_folder);


opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 1.5e3; %minimum atom number
opts.halo_N_lim = 10; %minimum allowed number in halo

opts.plot_dist = true; %do you want to see all the detailed stuff about the halo distributions

% % Background stuff
cli_header('Setting up for %s', data_folder);
opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
opts.data_src = fullfile(opts.data_root, data_folder);
opts.data_dir = data_folder;
opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
opts.index.filename = sprintf('index__%s__%.0f', opts.data_dir);
opts.label = data_folder;
opts.tag = 0;
opts.full_out = false;
opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds
opts.shot_bounds = [];
if ~exist(opts.fig_dir, 'dir')
    mkdir(opts.fig_dir);
end
% Run the function!

%% Set up out dir
%set up an output dir %https://gist.github.com/ferryzhou/2269380
if (exist([opts.data_src, '\out'], 'dir') == 0), mkdir(fullfile(opts.data_src, '\out')); end
%make a subfolder with the ISO timestamp for that date
anal_out.dir = sprintf('%sout\\%s\\', ...
    [opts.data_src, '\'], datestr(datetime('now'), 'yyyymmddTHHMMSS'));
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end

% % import raw data
[data, ~] = import_mcp_tdc_data(opts.import);
%% remove any ringing
opts.ring_lim = 0.5e-6; %how close can points be in time
data_masked = ring_removal(data,opts.ring_lim);
%% add labview import
if opts.tag
    logs = readtable(opts.logfile);
    tags = logs{:,5};
    %% select a specific shot type if you wish
    shot_type = 'double_halo';
    tag_mask = cellfun(@(x) strcmp(x, shot_type), tags');
    tag_mask = [tag_mask,zeros(1,length(data_masked.num_counts)-length(tags))];
else
    tag_mask = ones(1,length(data_masked.num_counts));
end
%% set up relevant constants
hebec_constants
%% find centers
opts.cent.visual = 0;
opts.cent.bin_size = 0.3e-5 * [1, 1, 1];
opts.cent.threshold = 2.0;
opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
bec = halo_cent(data_masked,opts.cent);
%% run some checks
% atoms number
% laser maybe?
num_check = data_masked.num_counts>opts.num_lim;
is_shot_good = num_check & bec.centre_OK_top' & bec.centre_OK_mid' & bec.centre_OK_btm' & tag_mask;
data_masked = struct_mask(data_masked,is_shot_good);
bec_masked = struct_mask(bec,is_shot_good);
%% Find the velocity widths
opts.bec_width.g0 = const.g0;
opts.bec_width.fall_time = 0.417;
bec_masked = bec_width_txy_to_vel(bec_masked,opts.bec_width);
%% convert data to velocity
%% generate top halo
opts.vel_conv.top.visual = 1;
opts.vel_conv.top.plot_percentage = 0.2;
opts.vel_conv.top.title = 'top halo';
opts.vel_conv.top.const.g0 = const.g0;
opts.vel_conv.top.const.fall_distance = const.fall_distance;
opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.top.v_mask=[0.9,1.09]; %bounds on radisu as multiple of radius value
opts.vel_conv.top.z_mask = [-0.04,0.04];

opts.vel_conv.top.bec_center.north = bec_masked.centre_top;
opts.vel_conv.top.bec_center.south = bec_masked.centre_mid;
opts.vel_conv.top.bec_width.north = bec_masked.width_top;
opts.vel_conv.top.bec_width.south = bec_masked.width_mid;
%%
top_halo = halo_vel_conv(data_masked,opts.vel_conv.top);
%% generate bottom halo
opts.vel_conv.btm.visual = 1;
opts.vel_conv.btm.plot_percentage = 0.2;
opts.vel_conv.btm.title = 'bottom halo';
opts.vel_conv.btm.const.g0 = const.g0;
opts.vel_conv.btm.const.fall_distance = const.fall_distance;
opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.btm.v_mask=[0.9,1.09]; %bounds on radisu as multiple of radius value
opts.vel_conv.btm.z_mask = [-0.04,0.04];

opts.vel_conv.btm.bec_center.north = bec_masked.centre_mid;
opts.vel_conv.btm.bec_center.south = bec_masked.centre_btm;
opts.vel_conv.btm.bec_width.north = bec_masked.width_mid;
opts.vel_conv.btm.bec_width.south = bec_masked.width_btm;
%%
bottom_halo = halo_vel_conv(data_masked,opts.vel_conv.btm);
%% mask out halos with nums to low
halo_N_check_top = top_halo.num_counts>opts.halo_N_lim;
halo_N_check_btm = bottom_halo.num_counts>opts.halo_N_lim;
halo_N_check = halo_N_check_top & halo_N_check_btm;
top_halo = struct_mask(top_halo,halo_N_check);
bottom_halo = struct_mask(bottom_halo,halo_N_check);
bec_masked = struct_mask(bec_masked,halo_N_check);
%% find the widths of the halos in velocity space
top_halo.bec_vel_width = (mean(bec_masked.vel_width_top,2)+mean(bec_masked.vel_width_mid,2))./2;% add the average bec width
bottom_halo.bec_vel_width = (mean(bec_masked.vel_width_btm,2)+mean(bec_masked.vel_width_mid,2))./2;
%% find the mode number
opts.mode_num.qe = 0.08;
top_halo.m = halo_mode_occupancy(top_halo,opts.mode_num);
bottom_halo.m = halo_mode_occupancy(bottom_halo,opts.mode_num);
%% calculated expected correlation amplitude
top_halo.g2 = 2 + 1./top_halo.m;
bottom_halo.g2 = 2 + 1./bottom_halo.m;
%% plot some histogram checks
halos.top_halo = top_halo;
halos.bottom_halo = bottom_halo;
halos.bec = bec_masked;
opts.plot_opts = [];
if opts.plot_dist
    plot_checks(halos,opts.plot_opts);
end
%% calculate correlation functions

%% back to back (intra halo)
corr_opts.fig='top halo bb corr';
corr_opts.type='1d_cart_bb';%'radial_bb';%'radial_cl';%'1d_cart_cl';%%
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*2e-3;
one_d_range=0.06;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,30)';
corr_opts.redges=(linspace(0,0.003,1500));
corr_opts.rad_smoothing=0.003;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;
corr_opts.plots=true;
corr_opts.norm_samp_factor=500;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;
corr_opts.calc_err=true;


% corr_opts.one_d_smoothing=nan;
corr_opts.one_d_smoothing=0.0008;

tic
top_halo.corr_bb=calc_any_g2_type(corr_opts,top_halo.counts_vel');
toc
%%
tic
corr_opts.fig='bottom halo bb corr';
bottom_halo.corr_bb=calc_any_g2_type(corr_opts,bottom_halo.counts_vel');
toc

% bb (inter halo)

%% co-linear (intra)
corr_opts.fig='top halo cl corr';
corr_opts.type='1d_cart_cl';%'radial_bb';%'radial_cl';%'1d_cart_cl';%%
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*10e-4;
one_d_range=0.03;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,50)';
corr_opts.redges=(linspace(0,0.003,1500));
corr_opts.rad_smoothing=0.003;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=nan;
corr_opts.plots=true;
corr_opts.norm_samp_factor=500;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;
corr_opts.calc_err=true;


corr_opts.one_d_smoothing=nan;
% corr_opts.one_d_smoothing=0.002;
tic
top_halo.corr_cl=calc_any_g2_type(corr_opts,top_halo.counts_vel');
toc
%%
tic
corr_opts.fig='bottom halo cl corr';
bottom_halo.corr_bb=calc_any_g2_type(corr_opts,bottom_halo.counts_vel');
toc

% cl (inter)