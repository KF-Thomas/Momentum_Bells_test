% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

path_config = 'C:\Users\jacob\Documents\Projects\QD\conf\config_20181213.m';
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% data_folder = '';
data_folder = '20191022_bell_prototyping';
opts.import.dir = fullfile(opts.data_root, data_folder);


opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

opts.num_lim = 1.5e3;

% % Background stuff
cli_header('Setting up for %s', data_folder);
opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
opts.data_src = fullfile(opts.data_root, data_folder);
opts.data_dir = data_folder;
opts = QD_config(opts);
opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
opts.index.filename = sprintf('index__%s__%.0f', opts.data_dir);
opts.label = data_folder;
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
%% set up relevant constants
hebec_constants
%% find centers
opts.cent.visual = 0;
opts.cent.bin_size = 0.5e-5 * [1, 1, 1];
opts.cent.threshold = 2.5;
t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972]}; %time bounds for the different momentum states k=+1,0,-1 respectively
%find the center of the top middle and bottom bec
opts.cent.crop = [t_bounds{1}; -0.03, 0.03; -0.03, 0.03];
[bec_centre_top,bec_width_top,bec_counts_top,centre_OK_top] =  find_dist_cen(data,opts.cent);
opts.cent.crop = [t_bounds{2}; -0.03, 0.03; -0.03, 0.03];
[bec_centre_mid,bec_width_mid,bec_counts_mid,centre_OK_mid] =  find_dist_cen(data,opts.cent);
opts.cent.crop = [t_bounds{3}; -0.03, 0.03; -0.03, 0.03];
[bec_centre_btm,bec_width_btm,bec_counts_btm,centre_OK_btm] =  find_dist_cen(data,opts.cent);
%% convert data to velocity
% top halo
for this_idx = 1:num_shots % Loop over all shots
    this_txy = data.counts_txy{this_idx};
    this_centre = (bec_centre_top(this_idx, :)+bec_centre_mid(this_idx, :))./2;
    centred_counts = this_txy - this_centre;

    % Convert to kspace
    this_outtime = this_centre(1) - 0.417;
    v_zxy = txy_to_vel(centred_counts, this_outtime, opts.const.g0, opts.const.fall_distance);
    top_halo.txy{this_idx} = this_txy;
    top_halo.vel{this_idx} = v_zxy;
end
% bottom halo