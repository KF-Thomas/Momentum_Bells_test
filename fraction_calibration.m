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
data_folder = '20191101_brief_halo_data';
opts.import.dir = fullfile(opts.data_root, data_folder);


opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;

opts.num_lim = 1.5e3;

% % Background stuff
cli_header('Setting up for %s', data_folder);
opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
opts.data_src = fullfile(opts.data_root, data_folder);
opts.data_dir = data_folder;
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
%% add labview import
% logs = readtable(opts.logfile);
%% run some checks
% atoms number
% laser maybe?
num_check = data.num_counts>opts.num_lim;
is_shot_good = num_check;
%% Count up in the different time bins
% time bounds we care about
t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972]}; %time bounds for the different momentum states k=+,1,0,-1 respectively
% t_bounds = {[3.8274,3.8611],[3.8611,3.9],[3.9,3.9467]}; %time bounds for the different magnetic states mj=+,1,0,-1 respectively

num_shots = length(data.shot_num);
Ns = zeros(num_shots,length(t_bounds));
for shot_idx = 1:num_shots
    for t_idx = 1:length(t_bounds)
        lims = [t_bounds{t_idx};opts.bounds];
        this_txy = data.counts_txy{shot_idx};
        trim_txy = masktxy_square(this_txy, lims);
        Ns(shot_idx,t_idx) = size(trim_txy,1);
    end
end
Ntotals = sum(Ns,2);
%% plot out the results
stfig('Fractions in each state')
plot(Ns./Ntotals)
ylabel('transfer fraction')
xlabel('Shot number')
legend('k=+1','k=0','k=-1')
set(gca,'fontsize',16)

