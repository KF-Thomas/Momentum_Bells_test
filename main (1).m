%% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
%% Import directory
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\kieran\Documents\LOCAL-DATA\';
% data_folder = '20191115_halos_attempt_3';
% data_folder = '20191114_halos_attempt_2';
% data_folder = '20200721_momentum_transfer_optimiser_2';
% data_folder = '20200721_momentum_transfer_optimiser_1';
% data_folder = '20200729_brief_movement_data_fat_cloud';
% data_folder = '20200807_k=0,-1,-2_halos_data_1';
% data_folder = '20200807_k=0,-1,-2_halos_data_2';
data_folder = '';
opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;
%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10

opts.plot_dist = true; %do you want to see all the detailed stuff about the halo distributions

%% Background stuff
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

%% import raw data
[data, ~] = import_mcp_tdc_data(opts.import);
%% remove any ringing
opts.ring_lim = 0.09e-6;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
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
opts.cent.visual = 2; %from 0 to 2
opts.cent.savefigs = 0;
opts.cent.correction = 1;
opts.cent.correction_opts.plots = 0;

opts.cent.top.visual = 0; %from 0 to 2
opts.cent.top.threshold = [130,2000,2000].*1e3;  %set in inverse units (Hz for time 1/m for space)
opts.cent.top.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.top.method = {'margin','average','average'};

opts.cent.mid.visual = 0; %from 0 to 2
opts.cent.mid.threshold = [130,2000,2000].*1e3;  %set in inverse units (Hz for time 1/m for space)
opts.cent.mid.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.mid.method = {'margin','average','average'};

opts.cent.btm.visual = 0; %from 0 to 2
opts.cent.btm.threshold = [130,2000,2000].*1e3;  %set in inverse units (Hz for time 1/m for space)
opts.cent.btm.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.btm.method = {'margin','average','average'};

% opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
% opts.cent.t_bounds = {[3.861,3.867],[3.874,3.881],[3.887,3.895],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]}; %time bounds for the different momentum states k=-2,-1,0 respectively
bec = halo_cent(data_masked,opts.cent);

%% run some checks
% atoms number
% laser maybe?
num_check = data_masked.num_counts>opts.num_lim;
% num_masked = data_masked.num_counts;
% num_masked(~num_check) = NaN;
% num_outlier = isoutlier(num_masked);
% ~num_outlier &
is_shot_good = num_check & bec.centre_OK_top' & bec.centre_OK_mid' & bec.centre_OK_btm' & tag_mask;
data_masked_halo = struct_mask(data_masked,is_shot_good);
bec_masked_halo = struct_mask(bec,is_shot_good);

%% Find the velocity widths
opts.bec_width.g0 = const.g0;
opts.bec_width.fall_time = 0.417;
bec_masked_halo = bec_width_txy_to_vel(bec_masked_halo,opts.bec_width);

%% convert data to velocity
% zero velocity point
t0 = ones(size(bec_masked_halo.centre_top,1),1).*3.8772;%bec_masked_halo.centre_top(:,1);%72;% 
x0 = bec_masked_halo.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
y0 = bec_masked_halo.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;

%% generate top halo
opts.vel_conv.top.visual = 0;
opts.vel_conv.top.plot_percentage = 0.95;
opts.vel_conv.top.title = 'top halo';
opts.vel_conv.top.const.g0 = const.g0;
opts.vel_conv.top.const.fall_distance = const.fall_distance;
opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.top.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
opts.vel_conv.top.z_mask = [-0.68,0.68]; %in units of radius (standard [-0.76,0.76])
opts.vel_conv.top.y_mask = [-1.9,1.9]; %in units of radius
opts.vel_conv.top.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point

opts.vel_conv.top.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s

opts.vel_conv.top.bec_center.north = bec_masked_halo.centre_top;
opts.vel_conv.top.bec_center.south = bec_masked_halo.centre_mid;
opts.vel_conv.top.bec_width.north = bec_masked_halo.width_top;
opts.vel_conv.top.bec_width.south = bec_masked_halo.width_mid;

%%
top_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.top);

%% generate bottom halo
opts.vel_conv.btm.visual = 0;
opts.vel_conv.btm.plot_percentage = 0.95;
opts.vel_conv.btm.title = 'bottom halo';
opts.vel_conv.btm.const.g0 = const.g0;
opts.vel_conv.btm.const.fall_distance = const.fall_distance;
opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.btm.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
opts.vel_conv.btm.z_mask = [-0.68,0.68]; %in units of radius
opts.vel_conv.btm.y_mask = [-1.9,1.9]; %in units of radius
opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point

opts.vel_conv.btm.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s

opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;

%%
bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);

%% mask out halos with nums to low
halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim;
halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim;
halo_N_check = halo_N_check_top & halo_N_check_btm;
top_halo = struct_mask(top_halo_intial,halo_N_check);
bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
bec_halo = struct_mask(bec_masked_halo,halo_N_check);

%% plot some histogram checks
halos.top_halo = top_halo;
halos.bottom_halo = bottom_halo;
halos.bec = bec_halo;
opts.plot_opts = [];
opts.plot_opts.only_dists = true;
opts.plot_opts.const.g0 = const.g0;
opts.plot_opts.const.fall_distance = const.fall_distance;
if opts.plot_dist
    plot_checks(halos,opts.plot_opts);
end

%% Separate data into the four ports
ports = {};
[ports.top_left, ports.top_right] = separate_ports(top_halo,0);
[ports.bottom_left, ports.bottom_right] = separate_ports(bottom_halo,0);

%% calculate the global correlation functions around the halos
global_corrs_opts.plots = true;
global_corrs_opts.fit = true;
global_corrs_opts.calc_err = false;

corrs = global_corrs(top_halo,bottom_halo,global_corrs_opts);

%% Quantum correlator E
opts_E.calc_err = false;
opts_E.plots = true;
opts_E.verbose = false;

[E_val, corrs.ports] = E(ports,opts_E);

%Expected amplitude
g14 = corrs.ports.g14.norm_g2.fitted_g2peak;
E_amp = (g14-1)/(g14+1);

%% Nice Plots of correlations
opts_nice_plots.calc_err = false;
opts_nice_plots.fit = true;

nice_g2_plots(opts_nice_plots,corrs.ports)

%% Write out results
top_halo_bb = string_value_with_unc(corrs.top_halo.corr_bb.norm_g2.fitted_g2peak,corrs.top_halo.corr_bb.norm_g2.fitted_g2peak_unc,'type','b');
top_halo_cl = string_value_with_unc(corrs.top_halo.corr_cl.norm_g2.fitted_g2peak,corrs.top_halo.corr_cl.norm_g2.fitted_g2peak_unc,'type','b');
bottom_halo_bb = string_value_with_unc(corrs.bottom_halo.corr_bb.norm_g2.fitted_g2peak,corrs.bottom_halo.corr_bb.norm_g2.fitted_g2peak_unc,'type','b');
bottom_halo_cl = string_value_with_unc(corrs.bottom_halo.corr_cl.norm_g2.fitted_g2peak,corrs.bottom_halo.corr_cl.norm_g2.fitted_g2peak_unc,'type','b');

between_halo_bb = string_value_with_unc(corrs.between_halos.corr_bb.norm_g2.fitted_g2peak,corrs.between_halos.corr_bb.norm_g2.fitted_g2peak_unc,'type','b');
between_halo_cl = string_value_with_unc(corrs.between_halos.corr_cl.norm_g2.fitted_g2peak,corrs.between_halos.corr_cl.norm_g2.fitted_g2peak_unc,'type','b');

cli_format_text('','c',3)
cli_format_text('RESULTS','c',3)
cli_format_text('','c',3)
fprintf('\n Top Halo:       g^{(2)}_{BB} = %s,    g^{(2)}_{CL} = %s\n',top_halo_bb,top_halo_cl)
fprintf('\n Bottom Halo:    g^{(2)}_{BB} = %s,    g^{(2)}_{CL} = %s\n',bottom_halo_bb,bottom_halo_cl)
fprintf('\n Between Halos:  g^{(2)}_{BB} = %s,    g^{(2)}_{CL} = %s\n',between_halo_bb,between_halo_cl)
fprintf('\n E = %g,   Expected Amplitude = %g\n',E_val,E_amp)


%HOM around the Halo
