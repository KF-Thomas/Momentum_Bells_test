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
% data_folder = '20191105_halos_3766_shots';
% data_folder = '20200623_mag_transfer_cal';
% data_folder = '20191104_halos_attempt_1';
% data_folder = '20191101_brief_halo_data';
% data_folder = '20200715_getting_bec_and_halo';
% data_folder = '20200717_mag_transfer_optimiser_2';
data_folder = '20200724_mag_transfer_optimiser_9';

frac_opts.num_lim = 0.4e3;%1.5e3;
frac_opts.transfer_state = 'momentum';%'mag';
frac_opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds 

% % Background stuff
cli_header('Setting up for %s', data_folder);
opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;
opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');

%% import raw data
[data, ~] = import_mcp_tdc_data(opts.import);
%% Count up in the different time bins
out_frac = fraction_calc(data,frac_opts);
%% plot out the results
stfig('Fractions')% in each momentum state: ',shot_type])
clf
scatter(out_frac.shot_num,out_frac.fracs(:,1))
hold on
scatter(out_frac.shot_num,out_frac.fracs(:,2))
scatter(out_frac.shot_num,out_frac.fracs(:,3))
ylabel('transfer fraction')
xlabel('Shot number')
title(['data: ',data_folder])
if strcmp(frac_opts.transfer_state,'momentum')
    legend('k=+1','k=0','k=-1')
elseif strcmp(frac_opts.transfer_state,'mag')
    legend('$m_J=+1$','$m_J=0$','$m_J=-1$')
end
set(gca,'fontsize',16)

