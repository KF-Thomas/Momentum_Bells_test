%% data vs phase

%% Initializing path
if ~exist('reload_data')
    clear all;
    loaded_data_folders = {};
end
    clear all;
% close all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
hebec_constants
combined_struct = @(S,T) cell2struct(cellfun(@vertcat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);

%% Import directories
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% opts.data_root = 'E:\he_bec_data_recovery\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\full_interferometer\rarity-tapster\tighter_trap\';
%data_folder = '';
log_folder = 'log_Phi.txt';
log_lab_folder = 'log_LabviewMatlab.txt';
data_folders = {
  %  '20230329_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_8'
     '20230325_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_5'
%   '20230324_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_extra'
%   '20230326_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_7'
%    '20230325_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_6'
 %    '20230324_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_4'
%     '20230324_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_3'
% '20230325_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S_5'

%     '20230323_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_S'

%     '20230323_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_4_improve'
%     '20230322_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_3_improve'
%     '20230322_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus_2'
%     '20230321_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_settings'

%     '20230317_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_amp_scan'
%     '20230317_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_amp_scan_1'
%     '20230317_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_amp_scan_2'
%     '20230318_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_amp_scan_3'
%     '20230319_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_amp_scan_4'

% '20230316_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus'
% '20230316_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_2'
% '20230317_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_200_mus_3'
% '20230320_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_170_mus'
% '20230321_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_250_mus'

%     '20230315_RT_k=+1,0,-1_Vsh_0_75_Vq_0_7_12000_mus'

% first go, only 2 points
%     '20211028_k=0,-1,-2_RT_run_1' %[pi/4,pi/2]
%     '20211028_k=0,-1,-2_RT_run_2' %[pi/4,pi/2]
%     '20211029_k=0,-1,-2_RT_run_3' %[pi/4,pi/2]
%     '20211029_k=0,-1,-2_RT_run_4' %[pi/4,pi/2]
    
%     '20211030_k=0,-1,-2_RT_run_5' %[0,pi/8,pi/4,pi/2]
%     '20211031_k=0,-1,-2_RT_run_6' %[0,pi/8,pi/4,pi/2]
%     '20211102_k=0,-1,-2_RT_run_7' %[0,pi/8,pi/4,3*pi/8,pi/2]
%     '20211103_k=0,-1,-2_RT_run_8' %[0,pi/8,pi/4,3*pi/8,pi/2]
%     '20211104_k=0,-1,-2_RT_run_9' %[0,pi/8,pi/4,3*pi/8,pi/2]
%     '20211105_k=0,-1,-2_RT_run_10'%[0,pi/8,pi/4,3*pi/8,pi/2]
%     '20211106_k=0,-1,-2_RT_run_11'%[0,pi/8,pi/4,3*pi/8,pi/2]
%     '20211107_k=0,-1,-2_RT_run_12'%[0,pi/8,pi/4,3*pi/8,pi/2]
    
% %     '20211108_k=0,-1,-2_RT_run_13'
%     %'20211109_k=0,-1,-2_RT_run_17' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2

%%%
%     '20211108_k=0,-1,-2_RT_run_14' %0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211109_k=0,-1,-2_RT_run_15' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211109_k=0,-1,-2_RT_run_16' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211109_k=0,-1,-2_RT_run_18'
%     '20211110_k=0,-1,-2_RT_run_19' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211110_k=0,-1,-2_RT_run_20' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211111_k=0,-1,-2_RT_run_21' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     
%     '20211117_adj_k=0,-1,-2_RT_run_23' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211118_adj_k=0,-1,-2_RT_run_24' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211118_adj_k=0,-1,-2_RT_run_25' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211119_adj_k=0,-1,-2_RT_run_26' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211119_adj_k=0,-1,-2_RT_run_27' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211120_adj_k=0,-1,-2_RT_run_28' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211122_adj_k=0,-1,-2_RT_run_29' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211123_adj_k=0,-1,-2_RT_run_30' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211124_adj_k=0,-1,-2_RT_run_31' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211126_adj_k=0,-1,-2_RT_run_32' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211129_adj_k=0,-1,-2_RT_run_33' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211130_adj_k=0,-1,-2_RT_run_34' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211130_adj_k=0,-1,-2_RT_run_35' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211130_adj_k=0,-1,-2_RT_run_36' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211130_adj_k=0,-1,-2_RT_run_37' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
%     '20211201_adj_k=0,-1,-2_RT_run_38'
%     '20211201_adj_k=0,-1,-2_RT_run_39'
%     '20211202_adj_k=0,-1,-2_RT_run_40'
%     '20211202_adj_k=0,-1,-2_RT_run_41'
%     '20211202_adj_k=0,-1,-2_RT_run_42'
%     '20211203_adj_k=0,-1,-2_RT_run_43'
%     '20211205_adj_k=0,-1,-2_RT_run_44'
%     '20211207_adj_k=0,-1,-2_RT_run_45'
%     '20211207_adj_k=0,-1,-2_RT_run_46'
%     '20211208_adj_k=0,-1,-2_RT_run_47'
%%%

    
    %first go
%     '20210803_k=0,-1,-2_rt_scan_mid_trap_equal_delay_1'
%     '20210804_k=0,-1,-2_rt_scan_mid_trap_equal_delay_2'%g34 not biased downward, good E
%     %     '20210805_k=0,-1,-2_rt_scan_mid_trap_equal_delay_3'% g34 phi=pi baised downward very strongly (not good E)
%     '20210806_k=0,-1,-2_rt_scan_mid_trap_equal_delay_4'% g34 phi=pi baised downward
%     '20210807_k=0,-1,-2_rt_scan_mid_trap_equal_delay_5'% g34 phi=pi baised downward but only slightly
%     '20210808_k=0,-1,-2_rt_scan_mid_trap_equal_delay_6'% g34 phi=pi baised downward
    

    %             '20210716_MZ_scan_50_mus_new_trap_14'%Max Vis = 0.54(4), at phi=0.115192
    %             '20210716_MZ_scan_50_mus_new_trap_13'%Max Vis = 0.58(5), at phi=0.0314159
    %             '20210716_MZ_scan_50_mus_new_trap_12'%Max Vis = 0.54(2), at phi=0.0314159
    %             '20210716_MZ_scan_50_mus_new_trap_11'%Max Vis = 0.46(3), at phi=0.115192
    %             '20210716_MZ_scan_50_mus_new_trap_10'%Max Vis = 0.43(2), at phi=-0.219911
    %             '20210713_MZ_scan_50_mus_new_trap_9' %Max Vis = 0.60(5), at phi=-0.219911
    %             '20210713_MZ_scan_50_mus_new_trap_8' %Max Vis = 0.55(5), at phi=-0.136136
    %             '20210713_MZ_scan_50_mus_new_trap_7' %Max Vis = 0.56(4), at phi=-0.136136
    %             '20210713_MZ_scan_50_mus_new_trap_6' %Max Vis = 0.58(4), at phi=-0.0523599
    %             '20210713_MZ_scan_50_mus_new_trap_5' %Max Vis = 0.52(3), at phi=-0.0523599
    %             '20210713_MZ_scan_50_mus_new_trap_4' %Max Vis = 0.59(4), at phi=-0.219911
    %             '20210713_MZ_scan_50_mus_new_trap_3' %Max Vis = 0.54(3), at phi=-0.219911
    %             '20210713_MZ_scan_50_mus_new_trap_2' %Max Vis = 0.592(10), at phi=-0.219911
    %             '20210713_MZ_scan_50_mus_new_trap'   %Max Vis = 0.57(3), at phi=-0.136136

    % RT attempts
    %     'full_interferometer\rarity-tapster\tighter_trap\20210525_RT_500_mus_length'
    %         'full_interferometer\rarity-tapster\20210219_rarity_tapster_k=0,-1,-2_scan_1200_mus_overnight'
    %         'full_interferometer\rarity-tapster\20210303_rarity_tapster_k=0,-1,-2_scan_1200_mus_overnight_v4'
    %         'full_interferometer\rarity-tapster\20210302_rarity_tapster_k=0,-1,-2_scan_1200_mus_overnight_v3'
    %          'full_interferometer\rarity-tapster\20210301_rarity_tapster_k=0,-1,-2_scan_1200_mus_overnight_v2'
    %         'full_interferometer\rarity-tapster\20210210_rarity_tapster_k=0,-1,-2_scan_2250_mus_overnight_v2'
    %         'full_interferometer\rarity-tapster\20210209_rarity_tapster_k=0,-1,-2_scan_2250_mus_overnight'
    %     '20210310_rarity_tapster_k=0,-1,-2_scan_1200_mus_evap_0_8543'
    %     '20210305_rarity_tapster_k=0,-1,-2_scan_1200_mus_evap_0_8543'
    
    
    %MZ attempts
    %     'full_interferometer\mach-zender\tighter_trap\20210604_MZ_scan_50_mus_new_pulses'
    %     'full_interferometer\mach-zender\tighter_trap\20210604_MZ_scan_50_mus_new_pulses_2'
    %     'full_interferometer\mach-zender\tighter_trap\20210604_MZ_scan_50_mus_new_pulses_3'
    %     'full_interferometer\mach-zender\tighter_trap\20210607_MZ_scan_50_mus_new_pulses_4'
    %     'full_interferometer\mach-zender\tighter_trap\20210607_MZ_scan_50_mus_new_pulses_5'
    %     'full_interferometer\mach-zender\20210224_mz_scan_1200_mus_shunt_0_65_attempts\20210224_mz_scan_1200_mus_length_chirped_pulses_shunt_0_65'
    %     'full_interferometer\mach-zender\20210224_mz_scan_1200_mus_shunt_0_65_attempts\20210224_mz_scan_1_2_ms_chirped_pulses_shunt_0_65_attempt_4'
    %     'full_interferometer\mach-zender\20210224_mz_scan_1200_mus_shunt_0_65_attempts\20210224_mz_scan_1_2_ms_chirped_pulses_shunt_0_65_attempt_2'
    %     'full_interferometer\mach-zender\20210225_mz_scan_1200_mus_length_shunt_0_65_tests\20210225_mz_scan_1200_mus_length_shunt_0_65_test_3'
    %     'full_interferometer\mach-zender\20210225_mz_scan_1200_mus_length_shunt_0_65_tests\20210225_mz_scan_1200_mus_length_shunt_0_65_test_4'
    %     'full_interferometer\mach-zender\20210219_mz_scan_1200_mus_length_chirped_pulses'
    % %     'full_interferometer\mach-zender\20210208_mach_zender_k=0,-1_distinguishability_dip\20210208_mach_zender_k=0,-1_scan_2250_mus'
    %         'full_interferometer\mach-zender\weak_trap\20210210_mach_zender_k=0,-1_lower_evap_long'
    %     'full_interferometer\mach-zender\20210219_mz_scan_1200_mus_length_chirped_pulses_attempt_2'
    %     'full_interferometer\mach-zender\20210218_mz_scan_1200_mus_length_opt_pulses'
    %     'full_interferometer\mach-zender\20210218_mz_scan_900_mus_length'
    %     'full_interferometer\mach-zender\20210223_mz_scan_1200_mus_length_chirped_pulses_after_coil_fail'
    %     '20210310_mz_scan_1300_mus_length'
    %     '20210315_mz_scan_k=-1,-2_1200_mus_length_improved_pulses'
    %     '20210315_mz_scan_k=0,-1_1200_mus_length_improved_pulses'
    };

force_reimport = false;
force_reimport_norm = true;

%%
if exist('data','var') && cell_comp(data_folders, loaded_data_folders)
   reload_data = false;
   clear out_data phi_vec out_corrs
else
   reload_data = true;
   loaded_data_folders = data_folders;
   clear data out_data phi_vec out_corrs
end

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 3.53e3;%1.5e3;%2.3e3;%9e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper = Inf;%6;%8;%15;%10;%17;%12;%30;%5;%10;%20;%3.5;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_norm = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper_norm = Inf;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_both = 1;%1;
opts.halo_number_limit = Inf;%24;%25;%Inf;%26;%30;%21;

y_cut = 11e-3;
catch_count = 0;

ang_lim = 30;%65;%23;%55;%13;%13;%11;%15;%20.0;%20;%40;%%20;%angular limit in degrees



% opts.import.shot_num = 1:16; %can select which shots you want to import

%% Calibration settings
L = [-0.3,0.3];%[0,0.3].*2;%[-0.4,-0.12].*2;%[-0.15,0.15];%[-0.1,0.1];%[-0.618,-0.385];%[-0.5,-0.1];%
L_theta =[-2.*pi, 2.*pi];%[0.5,1.0].*pi;%[-0.5,0.5];%
norm_folders = [];%[];%
plot_fit = true;
do_g2=true;
do_g2_err=false;
do_E_err=false;
do_range_cut=false;
out_data = {};
phi_vec = [];

num_shots_norm = 0;
opts.tag = 0;

%% Run over each folder
% close all
for folder_indx = 1:length(data_folders)
    %% import raw data
    data_folder = data_folders{folder_indx};
    opts.import.dir = fullfile(opts.data_root, data_folder);
    opts.logfile = fullfile(opts.import.dir,log_folder);
    opts.loglabfile = fullfile(opts.import.dir,log_lab_folder);
    if ~ismember(folder_indx,norm_folders)
        phi_logs = table2array(readtable(opts.logfile));
        opts.cent.nan_cull = true;
        opts.import.force_reimport = force_reimport;
        opts.import.force_cache_load = ~opts.import.force_reimport;
        %         opts.import.shot_num = 1:24; %can select which shots you want to import
    else
        opts.cent.nan_cull = false;
        opts.import.force_reimport = force_reimport_norm;
        opts.import.force_cache_load = ~opts.import.force_reimport;
    end
    
    
    %% Chose which halo(s) to analyse
    opts.do_top_halo = 1;% analyse the top halo?
    opts.do_btm_halo = 1;% analyse the bottom halo?
    
    %% Chose if you want to look at a narrow or wide slice of the halo
    slice_type = 'extra narrow';
    if strcmp(slice_type,'narrow')
        opts.vel_conv.top.z_mask = [-0.4,0.4];%
        opts.vel_conv.btm.z_mask = [-0.4,0.4];%in units of radius ([-0.68,0.68])
    elseif strcmp(slice_type,'medium')
        opts.vel_conv.top.z_mask = [-0.6,0.6];%
        opts.vel_conv.btm.z_mask = [-0.6,0.6];%in units of radius ([-0.68,0.68])
    elseif strcmp(slice_type,'extra wide')
        opts.vel_conv.top.z_mask = [-0.87,0.87];%
        opts.vel_conv.btm.z_mask = [-0.87,0.87];%in units of radius ([-0.68,0.68])
    elseif strcmp(slice_type,'extra narrow')
        opts.vel_conv.top.z_mask = [-0.2,0.2];%
        opts.vel_conv.btm.z_mask = [-0.2,0.2];%in units of radius ([-0.68,0.68])
    elseif strcmp(slice_type,'super narrow')
        opts.vel_conv.top.z_mask = [-0.1,0.1];%
        opts.vel_conv.btm.z_mask = [-0.1,0.1];%in units of radius ([-0.68,0.68])
    else
        opts.vel_conv.top.z_mask = [-0.82,0.82];
        opts.vel_conv.btm.z_mask = [-0.82,0.82];%in units of radius ([-0.68,0.68])
    end
    
    opts.vel_conv.top.z_mask = [-0.9,0.9];
    opts.vel_conv.btm.z_mask = [-0.9,0.9];%in units of radius ([-0.68,0.68])
    radius_lim = [0.057,0.075];%[0.0515,0.079];%[0.062,0.071];%[0.03,0.09];%[0.03,0.09];%[0.06,0.069];%[0.03,0.09];%[0.0585,0.069];%[0.06,0.067];%[0.06,0.07];%[0.045,0.085];%[0.79,1.17].*0.065;%[0.61,1.26];%[0.89,1.11];%[0.89,1.16];%[0.9,1.05];%
    %     ang_lim = 20;%40;%%20;%angular limit in degrees
    %     ang_lim = 90;%
    
    %% Import parameters
    tmp_xlim=[-55e-3, 55e-3]; %[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
    tmp_ylim=[-55e-3, 55e-3]; %[-35e-3, 35e-3];
    tlim=[0,4];
    opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];
    
%     opts.num_lim = num_lim;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
    
    opts.plot_dist = false; %do you want to see all the detailed stuff about the halo distributions
    
    %% Background stuff
    cli_header('Setting up for %s', data_folder);
    opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
    opts.data_src = fullfile(opts.data_root, data_folder);
    opts.data_dir = data_folder;
    opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
    %     opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
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
    
    if reload_data
    %% import raw data
    [data{folder_indx}, ~] = import_mcp_tdc_data(opts.import);
    
    %% remove any ringing or hotspots
%     data_ht_spot=hotspot_mask(data{folder_indx});
%     data{folder_indx}.counts_txy=data_ht_spot.masked.counts_txy;
%     data{folder_indx}.num_counts=data_ht_spot.masked.num_counts;
    end
    
    opts.ring_lim = 0.09e-6;%-1;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
    data_masked = ring_removal(data{folder_indx},opts.ring_lim);
    if folder_indx == 1
       stfig('num history');
       clf
       hold on
       plot(data_masked.num_counts)
    else
        stfig('num history');
        plot(data_masked.num_counts)
    end
    
    %% add labview import
    if opts.tag && ~ismember(folder_indx,norm_folders)
        opts_tab = detectImportOptions(opts.loglabfile);
        opts_tab.Delimiter = {','};
        logs = readtable(opts.loglabfile,opts_tab);
        tags = logs{:,6};
        %% select a specific shot type if you wish
        shot_type = 'main';
        tag_mask = cellfun(@(x) strcmp(x, shot_type), tags');
        tag_mask = [tag_mask,zeros(1,length(data_masked.num_counts)-length(tags))];
        data_masked = struct_mask(data_masked,logical(tag_mask),1);
        logs_cell = readcell(opts.loglabfile,opts_tab);
        
        settings_log=cell2mat(logs_cell(:,5));
        settings_log=settings_log(:,18:end);
        settings_list = settings_list_func();
        [Tf, Loc]=ismember(settings_log, settings_list);
        evap_setting = Loc(logical(tag_mask));
    end
    
    
    %% set up relevant constants
    hebec_constants
    
    %% find centers
    opts.cent.visual = 0; %from 0 to 2
    opts.cent.savefigs = 0;
    opts.cent.correction = 0;%1;
    opts.cent.correction_opts.plots = 0;
    
    opts.cent.top.visual = 0; %from 0 to 2
    opts.cent.top.savefigs = 0;
    opts.cent.top.threshold = [150,6000,6000].*1e3;%130 vs 70 % [130,6000,6000].*1e3
    opts.cent.top.min_threshold = [0,20,20].*1e3;%[0,3,3].*1e3;%[16,7,10].*1e3; [0,30,30]  [0,8,8]
    opts.cent.top.sigma = [6.7e-5,16e-5,17e-5];%[6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.top.method = {'margin','average','average'};%{'margin','average','average'};
    
    opts.cent.mid.visual = 0; %from 0 to 2
    opts.cent.mid.savefigs = 0;
    opts.cent.mid.threshold = [150,6000,6000].*1e3;%130
    opts.cent.mid.min_threshold = [0,20,20].*1e3;%[16,7,10].*1e3;
    opts.cent.mid.sigma = [6.7e-5,16e-5,17e-5];%[8e-5,25e-5,25e-5];
    opts.cent.mid.method = {'margin','average','average'};
    
    opts.cent.btm.visual = 2; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [150,6000,6000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [0,20,20].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.btm.sigma = [6.7e-5,16e-5,17e-5];%[8e-5,25e-5,25e-5];
    opts.cent.btm.method = {'margin','average','average'};
    
    %              opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
        opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.884,3.896],[3.75,4]};%

%     opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
    bec = halo_cent(data_masked,opts.cent);
    
    %%
    l_t=phi_logs(:,2);
    %     d_t=data_masked.time_create_write(:,1);
    d_t=data_masked.time_create_write(:,2);
    phi_log_matched = zeros(length(d_t),1).*nan;
    phi_log_check = zeros(length(d_t),1);
    for ii = 1: length(l_t)
        phi_c = phi_logs(ii,3);
        l_c=l_t(ii);
        t_mask=l_c+17<d_t & l_c+31>d_t;%29.5
        d_indx=find(t_mask);
        phi_log_matched(d_indx,1) = phi_c;
        phi_log_check(d_indx,1) = 1;
        
    end
    
    %% run some checks
    % atoms number
    % laser maybe?
    num_check = data_masked.num_counts>opts.num_lim;
    % num_masked = data_masked.num_counts;
    % num_masked(~num_check) = NaN;
    % num_outlier = isoutlier(num_masked);
    % ~num_outlier &
    is_shot_good = num_check & bec.centre_OK_mid';% & phi_log_check';
    if ~ismember(folder_indx,norm_folders)
        is_shot_good = is_shot_good & phi_log_check';
    end
    if ~ismember(folder_indx,norm_folders) && do_range_cut
        slosh_cut = ~(abs(bec.centre_mid(:,3))>y_cut);
        is_shot_good = slosh_cut' & is_shot_good;
    end
    if opts.do_top_halo
        is_shot_good = is_shot_good & bec.centre_OK_top';
    end
    if opts.do_btm_halo
        is_shot_good = is_shot_good & bec.centre_OK_btm';
    end
    data_masked_halo = struct_mask(data_masked,is_shot_good);
    bec_masked_halo = struct_mask(bec,is_shot_good);
    if opts.tag
        evap_setting = evap_setting(is_shot_good);
    end
    if ~ismember(folder_indx,norm_folders)
        phi_logs_masked = phi_log_matched(is_shot_good,:);
    end
    %% Find the velocity widths
    opts.bec_width.g0 = const.g0;
    opts.bec_width.fall_time = 0.417;
    bec_masked_halo = bec_width_txy_to_vel(bec_masked_halo,opts.bec_width);
    
    %% convert data to velocity
    % zero velocity point
    t0 = bec_masked_halo.centre_top(:,1);%ones(size(bec_masked_halo.centre_top,1),1).*3.8772;%72;%
    x0 = bec_masked_halo.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
    y0 = bec_masked_halo.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;
    
    %% generate top halo
    opts.vel_conv.top.visual = 0;
    opts.vel_conv.top.plot_percentage = 0.95;
    opts.vel_conv.top.title = 'top halo';
    opts.vel_conv.top.const.g0 = const.g0;
    opts.vel_conv.top.const.fall_distance = const.fall_distance;
    opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
    opts.vel_conv.top.v_mask=radius_lim; %bounds on radisu as multiple of radius value
    opts.vel_conv.top.ang_lim = ang_lim; %angular limits of the azimuthal angle
    opts.vel_conv.top.y_mask = [-1.9,1.9]; %in units of radius
    opts.vel_conv.top.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point
    
    opts.vel_conv.top.theta_mask = [];%[1.4,2.6;1.4-pi,2.6-pi];%[];%[0.6 2.45;0.6-pi 2.45-pi];%[1.17,1.28;1.17-pi,1.28-pi];%[];%[1.17,1.28;1.17-pi,1.28-pi];%

    opts.vel_conv.top.theta_mask_polarity = 0;

    opts.vel_conv.top.centering_correction = [0,0,0]; %[1.769005720082404e+00     1.058658203190879e-01    -3.319223909672235e-01].*0.5e-3;%correctoin shift to the centering in m/s
    
    opts.vel_conv.top.bec_center.north = bec_masked_halo.centre_top;
    opts.vel_conv.top.bec_center.south = bec_masked_halo.centre_mid;
    opts.vel_conv.top.bec_width.north = bec_masked_halo.width_top;
    opts.vel_conv.top.bec_width.south = bec_masked_halo.width_mid;
    
    %%
    if opts.do_top_halo
        top_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.top);
        C    = cell(1, size(top_halo_intial.counts_txy,2));
        C(:) = {data_folder};
        top_halo_intial.data_folder = C;
    else
        top_halo_intial.counts_vel = {};
        top_halo_intial.counts_vel_norm = {};
        top_halo_intial.num_counts = [];
    end
    
    %% generate bottom halo
    opts.vel_conv.btm.visual = 0;
    opts.vel_conv.btm.plot_percentage = 0.95;
    opts.vel_conv.btm.title = 'bottom halo';
    opts.vel_conv.btm.const.g0 = const.g0;
    opts.vel_conv.btm.const.fall_distance = const.fall_distance;
    opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
    opts.vel_conv.btm.v_mask=radius_lim; %bounds on radisu as multiple of radius value
    opts.vel_conv.btm.ang_lim = ang_lim; %angular limits of the azimuthal angle
    opts.vel_conv.btm.y_mask = [-1.9,1.9]; %in units of radius
    opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point
    
    opts.vel_conv.btm.theta_mask = [];%[1.4,2.6;1.4-pi,2.6-pi];%[0.6 2.45;0.6-pi 2.45-pi];%;%[1.17,1.28;1.17-pi,1.28-pi];%[];%[1.17,1.28;1.17-pi,1.28-pi];%

    opts.vel_conv.btm.theta_mask_polarity = 0;

    opts.vel_conv.btm.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s
    opts.vel_conv.btm.phi_correction = [0,0];

    opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
    opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
    opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
    opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;
    
    %%
    if opts.do_btm_halo
        bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);
        C    = cell(1, size(bottom_halo_intial.counts_txy,2));
        C(:) = {data_folder};
        bottom_halo_intial.data_folder = C;
    else
        bottom_halo_intial.counts_vel = {};
        bottom_halo_intial.counts_vel_norm = {};
        bottom_halo_intial.num_counts = [];
    end
    
    %% mask out halos with nums to low
    if ismember(folder_indx,norm_folders)
        halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim_norm & top_halo_intial.num_counts<opts.halo_N_lim_upper_norm;
        halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim_norm & bottom_halo_intial.num_counts<opts.halo_N_lim_upper_norm;
        halo_N_check_both = halo_N_check_top & halo_N_check_btm;
        
        full_halo_N_check_both= halo_N_check_both;
    else
        halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim & top_halo_intial.num_counts<opts.halo_N_lim_upper;
        halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim & bottom_halo_intial.num_counts<opts.halo_N_lim_upper;
        halo_N_check_both = (top_halo_intial.num_counts+bottom_halo_intial.num_counts)>opts.halo_N_lim_both;
        
        full_halo_N_check_both = (top_halo_intial.top_halo_num+bottom_halo_intial.btm_halo_num)<opts.halo_number_limit;
        
    end
    if opts.do_top_halo && opts.do_btm_halo
        halo_N_check = halo_N_check_top & halo_N_check_btm &halo_N_check_both&full_halo_N_check_both;
        top_halo = struct_mask(top_halo_intial,halo_N_check);
        bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
    elseif opts.do_top_halo
        halo_N_check = halo_N_check_top;
        top_halo = struct_mask(top_halo_intial,halo_N_check);
        bottom_halo = bottom_halo_intial;
    else
        halo_N_check = halo_N_check_btm;
        top_halo = top_halo_intial;
        bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
    end
    bec_halo = struct_mask(bec_masked_halo,halo_N_check);
    if opts.tag
    evap_setting = evap_setting(halo_N_check);
    end
    
    %% Seperate halos by phase
    if ismember(folder_indx,norm_folders)
        %% histograming
        nbins=151;%50;%
        theta_bins = linspace(-pi,pi,nbins+1);
        phi_bins = linspace(-pi/2,pi/2,nbins+1);
        v_top_zxy = cell2mat(top_halo.counts_vel_norm);
        v_top_zxy_unnorm = cell2mat(top_halo.counts_vel);
        v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm);
        v_btm_zxy_unnorm = cell2mat(bottom_halo.counts_vel);
        if opts.do_top_halo
            [theta_top,rxy_top] = cart2pol(v_top_zxy_unnorm(:,2),v_top_zxy_unnorm(:,3));
            phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
        else
            theta_top = [];
            phi_top = [];
        end
        
        if opts.do_btm_halo
            [theta_btm,rxy_btm] = cart2pol(v_btm_zxy_unnorm(:,2),v_btm_zxy_unnorm(:,3));
            phi_btm = atan(v_btm_zxy(:,1)./sqrt(v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2));
        else
            theta_btm = [];
            phi_btm = [];
        end
        
        v_btm_dens = [];
        v_top_dens = [];
        v_btm_dens_unc = [];
        v_top_dens_unc = [];
        %     phi_mask_top = (phi_top<0.154& phi_top>-0.154);
        %     phi_mask_btm = (phi_btm<0.154& phi_btm>-0.154);
        %     r_btm_zxy_masked=smooth_hist(theta_btm(phi_mask_btm),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
        %     r_top_zxy_masked=smooth_hist(theta_top(phi_mask_top),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
        %     v_btm_dens(:,1) = r_btm_zxy_masked.count_rate.smooth;
        %     v_top_dens(:,1) = r_top_zxy_masked.count_rate.smooth;
        %     v_btm_dens_unc(:,1) = sqrt(r_btm_zxy_masked.count_rate.smooth).*sqrt(abs(r_btm_zxy_masked.bin.edge(1:end-1)...
        %         -r_btm_zxy_masked.bin.edge(2:end)));
        %     v_top_dens_unc(:,1) = sqrt(r_top_zxy_masked.count_rate.smooth).*sqrt(abs(r_top_zxy_masked.bin.edge(1:end-1)...
        %         -r_top_zxy_masked.bin.edge(2:end)));
        
        
        r_btm_zxy_masked=smooth_hist(phi_btm,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
        r_top_zxy_masked=smooth_hist(phi_top,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
        v_btm_dens(:,2) = r_btm_zxy_masked.count_rate.smooth;
        v_top_dens(:,2) = r_top_zxy_masked.count_rate.smooth;
        
        if ~exist('cal_dens_top','var')
            cal_dens_top = zeros(nbins,1);
            cal_dens_btm = zeros(nbins,1);
        end
        cal_dens_top = cal_dens_top + v_top_dens(:,2)./size(norm_folders,1);
        cal_dens_btm = cal_dens_btm + v_btm_dens(:,2)./size(norm_folders,1);
        num_shots_norm = num_shots_norm+size(top_halo.counts_vel_norm,1);
    else
        
        unique_phi = unique(phi_logs_masked(halo_N_check,1)); %unique phases used in this scan
        unique_phi = unique_phi(~isnan(unique_phi));
        for ii = 1:length(unique_phi)
            current_phi = unique_phi(ii,1);
            phi_mask = (phi_logs_masked(halo_N_check,1)==current_phi);
            masked_top = struct_mask(top_halo,phi_mask');
            masked_bottom = struct_mask(bottom_halo,phi_mask');
            if opts.tag
                masked_evap = evap_setting(phi_mask);
            end
            if ~ismember(current_phi,phi_vec)
                phi_vec = [phi_vec,current_phi];
            end
            phi_indx = find(phi_vec==current_phi);
            if length(out_data)<phi_indx
                out_data{phi_indx} = {};
                out_data{phi_indx}.top_halo = masked_top;
                out_data{phi_indx}.bottom_halo = masked_bottom;
%                 out_data{phi_indx}.evap_setting = masked_evap;
            else
                try
                out_data{phi_indx}.top_halo = combined_struct(out_data{phi_indx}.top_halo,masked_top);
                out_data{phi_indx}.bottom_halo =combined_struct(out_data{phi_indx}.bottom_halo,masked_bottom);
%                 out_data{phi_indx}.evap_setting = masked_evap;
                catch
                    catch_count = catch_count+1;
                end
            end
        end
    end
    %% looking at sloshing
    %     nanstd(bec.centre_top(:,1))
    %     nanstd(bec.centre_top(:,2))
    %     nanstd(bec.centre_top(:,3))
    
end
%
%% Setting up variables
%     halos.top_halo = top_halo;
%     halos.bottom_halo = bottom_halo;
%     halos.bec = bec_halo;
opts.plot_opts = [];
opts.plot_opts.only_dists = true;
opts.plot_opts.const.g0 = const.g0;
opts.plot_opts.const.fall_distance = const.fall_distance;
% if opts.plot_dist
%     plot_checks(halos,opts.plot_opts);
% end
top_dens_vec = [];
top_dens_ind = [];
out_conv_dens = zeros(length(phi_vec),76);

if exist('cal_dens_top','var')
    cal_dens_top = cal_dens_top./num_shots_norm;
    cal_dens_btm = cal_dens_btm./num_shots_norm;
end
E_lambda = [];

for ii = 1:length(phi_vec)
    d = opts.plot_opts.const.fall_distance;
    g0 = opts.plot_opts.const.g0;
    tf = sqrt(2*d/g0);%arrival time of zero velocity particles
    
    v_top_zxy = cell2mat(out_data{ii}.top_halo.counts_vel_norm);
    v_top_zxy_unnorm = cell2mat(out_data{ii}.top_halo.counts_vel);
    r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
    r_dist_top_unnorm = sqrt(v_top_zxy_unnorm(:,1).^2+v_top_zxy_unnorm(:,2).^2+v_top_zxy_unnorm(:,3).^2);
    N_top = top_halo.num_counts;
    
    v_btm_zxy = cell2mat(out_data{ii}.bottom_halo.counts_vel_norm);
    v_btm_zxy_unnorm = cell2mat(out_data{ii}.bottom_halo.counts_vel);
    if ~isempty(v_btm_zxy)
        r_dist_btm = sqrt(v_btm_zxy(:,1).^2+v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2);
        r_dist_btm_unnorm = sqrt(v_btm_zxy_unnorm(:,1).^2+v_btm_zxy_unnorm(:,2).^2+v_btm_zxy_unnorm(:,3).^2);
    else
        r_dist_btm = [];
        r_dist_btm_unnorm = [];
    end
    N_btm = bottom_halo.num_counts;
    
    %% histograming
    nbins=151;%50;%
    theta_bins = linspace(-pi,pi,nbins+1);
    phi_bins = linspace(-pi/2,pi/2,nbins+1);
    if opts.do_top_halo
        [theta_top,rxy_top] = cart2pol(v_top_zxy_unnorm(:,2),v_top_zxy_unnorm(:,3));
        phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
    else
        theta_top = [];
        phi_top = [];
    end
    
    if opts.do_btm_halo
        [theta_btm,rxy_btm] = cart2pol(v_btm_zxy_unnorm(:,2),v_btm_zxy_unnorm(:,3));
        phi_btm = atan(v_btm_zxy(:,1)./sqrt(v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2));
    else
        theta_btm = [];
        phi_btm = [];
    end
    num_shots = size(out_data{ii}.top_halo.counts_vel_norm,1);
    v_btm_dens = [];
    v_top_dens = [];
    v_btm_dens_unc = [];
    v_top_dens_unc = [];
    phi_mask_top = (phi_top<0.154& phi_top>-0.154);
    phi_mask_btm = (phi_btm<0.154& phi_btm>-0.154);
    r_btm_zxy_masked=smooth_hist([theta_btm(phi_mask_btm)-2*pi;theta_btm(phi_mask_btm);theta_btm(phi_mask_btm)+2*pi],'sigma',0.12,'lims',[-pi,pi],'bin_num',nbins);
    r_top_zxy_masked=smooth_hist([theta_top(phi_mask_top)-2*pi;theta_top(phi_mask_top);theta_top(phi_mask_top)+2*pi],'sigma',0.12,'lims',[-pi,pi],'bin_num',nbins);
    theta_bins = r_btm_zxy_masked.bin.centers;
    v_btm_dens(:,1) = r_btm_zxy_masked.count_rate.smooth./num_shots;
    v_top_dens(:,1) = r_top_zxy_masked.count_rate.smooth./num_shots;
    v_btm_dens_unc(:,1) = sqrt(r_btm_zxy_masked.count_rate.smooth).*sqrt(abs(r_btm_zxy_masked.bin.edge(1:end-1)...
        -r_btm_zxy_masked.bin.edge(2:end)));
    v_top_dens_unc(:,1) = sqrt(r_top_zxy_masked.count_rate.smooth).*sqrt(abs(r_top_zxy_masked.bin.edge(1:end-1)...
        -r_top_zxy_masked.bin.edge(2:end)));
    
    theta_mask_btm = theta_btm<L_theta(2) & theta_btm>L_theta(1);
    theta_mask_top = theta_top<L_theta(2) & theta_top>L_theta(1);
    r_btm_zxy_masked=smooth_hist(phi_btm(theta_mask_btm),'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
    r_top_zxy_masked=smooth_hist(phi_top(theta_mask_top),'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
    phi_bins = r_btm_zxy_masked.bin.centers;
    v_btm_dens(:,2) = r_btm_zxy_masked.count_rate.smooth./num_shots;
    v_top_dens(:,2) = r_top_zxy_masked.count_rate.smooth./num_shots;
    
    v_top_dens_2d = hist3([theta_top phi_top],'Nbins',[nbins nbins]);
    v_btm_dens_2d = hist3([theta_btm phi_btm],'Nbins',[nbins nbins]);
    
    theta = linspace(-pi,pi,nbins);
    phi = linspace(-pi/2,pi/2,nbins);
    %%
    stfig('radial distribution phase');
    % clf
    r_hist_top=smooth_hist(r_dist_top,'sigma',0.0001,'lims',[0.03,0.07],'bin_num',600);
    r_hist_btm=smooth_hist(r_dist_btm,'sigma',0.0001,'lims',[0.03,0.07],'bin_num',600);
    r_hist_top_un=smooth_hist(r_dist_top_unnorm,'sigma',0.0001,'lims',[0.03,0.08],'bin_num',600);
    r_hist_btm_un=smooth_hist(r_dist_btm_unnorm,'sigma',0.0001,'lims',[0.03,0.08],'bin_num',600);
    subplot(2,1,1)
    plot(r_hist_top_un.bin.centers,r_hist_top_un.counts.smooth./(r_hist_top_un.counts.smooth+r_hist_btm_un.counts.smooth),'linewidth',1.5)
    hold on
    % plot(r_hist_top.bin.centers,r_hist_top.counts.smooth,'linewidth',1.5)
    % hold on
    % plot(r_hist_btm.bin.centers,r_hist_btm.counts.smooth,'linewidth',1.5)
    % xlabel('r')
    % ylabel('Freq')
    % xlim([min([r_hist_top.bin.centers;r_hist_btm.bin.centers]),...
    %     max([r_hist_top.bin.centers;r_hist_btm.bin.centers])])
    % legend('top','btm')
    subplot(2,1,2)
    plot(r_hist_top_un.bin.centers,r_hist_top_un.counts.smooth,'linewidth',1.5)
    hold on
    plot(r_hist_btm_un.bin.centers,r_hist_btm_un.counts.smooth,'linewidth',1.5)
    ylimit = max([r_hist_top_un.counts.smooth;r_hist_btm_un.counts.smooth]);
    % plot([0.130159/2 0.130159/2],[-0.1,ylimit.*2],'k-','linewidth',1.5)
    % ylim([0 ylimit.*1.1])
    xlabel('r')
    ylabel('Freq')
    xlim([min([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers]),...
        max([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers])])
    legend('top','btm','expected radius')
    
    %%
    for jj = 1:size(out_data{ii}.top_halo.counts_vel,1)
        current_counts = out_data{ii}.top_halo.counts_vel{jj};
        curret_num_counts = out_data{ii}.top_halo.num_counts(jj) + out_data{ii}.bottom_halo.num_counts(jj);
        current_shot_num = out_data{ii}.top_halo.shot_num(jj);
        current_counts_btm = out_data{ii}.bottom_halo.counts_vel{jj};
        phi_c = atan(current_counts(:,1)./sqrt(current_counts(:,2).^2+current_counts(:,3).^2));
        phi_c_btm = atan(current_counts_btm(:,1)./sqrt(current_counts_btm(:,2).^2+current_counts_btm(:,3).^2));
        [theta_c,rxy_c] = cart2pol(current_counts(:,2),current_counts(:,3));
        phi_mask_c = (phi_c<L(2)/2& phi_c>L(1)/2);
        if isempty(theta_c(phi_mask_c))
            continue
        end
        theta_dens_c = smooth_hist(theta_c(phi_mask_c),'sigma',0.1,'lims',[-pi,pi],'bin_num',nbins);
        
        
        pi_indx = floor(nbins/2);
        conv_theta_dens_c = theta_dens_c.count_rate.smooth(1:end-pi_indx).*theta_dens_c.count_rate.smooth(pi_indx+1:end);
        out_conv_dens(ii,:) = out_conv_dens(ii,:) + conv_theta_dens_c'./size(out_data{ii}.top_halo.counts_vel,1);
        
        r_btm_zxy_masked=smooth_hist(phi_c_btm,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
        r_top_zxy_masked=smooth_hist(phi_c,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
        v_btm_dens_c = r_btm_zxy_masked.count_rate.smooth;
        v_top_dens_c = r_top_zxy_masked.count_rate.smooth;
        
        %% ratio in spherical coordinates
        top_dens_norm_c = (v_top_dens_c)./(v_btm_dens_c+v_top_dens_c);
        phi_mask = phi<L(2)/2 & phi>L(1)/2;
        top_dens_avg_c = trapz(phi(phi_mask),top_dens_norm_c(phi_mask))./range(phi(phi_mask));
        top_dens_std_c = sqrt(trapz(phi(phi_mask),(top_dens_avg_c-top_dens_norm_c(phi_mask)).^2)./(range(phi(phi_mask))));
        
        top_dens_ind = [top_dens_ind;[phi_vec(ii),top_dens_avg_c,top_dens_std_c,curret_num_counts,current_shot_num]];
    end
    %     cell2mat(out_data{ii}.top_halo.counts_vel);
    %     smooth_hist(theta_btm(phi_mask_btm),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
    
    
    %% ratio in spherical coordinates
    top_dens_norm = (v_top_dens(:,:))./(v_btm_dens(:,:)+v_top_dens(:,:));
    phi_mask = phi_bins<L(2)/2 & phi_bins>L(1)/2;
    theta_mask = theta_bins<L_theta(2) & theta_bins>L_theta(1);
    ang_mask = phi_mask;% & theta_mask;
    top_dens_avg = trapz(phi(ang_mask),top_dens_norm(ang_mask,2))./range(phi(ang_mask));
    top_dens_std = sqrt(trapz(phi(ang_mask),(top_dens_avg-top_dens_norm(ang_mask,2)).^2)./(range(phi(phi_mask))));
    %
    %     top_dens_avg = trapz(theta(theta_mask),top_dens_norm(theta_mask,2))./range(theta(theta_mask));
    %     top_dens_std = sqrt(trapz(theta(theta_mask),(top_dens_avg-top_dens_norm(theta_mask,2)).^2)./(range(theta(theta_mask))));
    
    if ~exist('cal_dens_top','var')
        cal_dens_top = zeros(nbins,1);
        cal_dens_btm = zeros(nbins,1);
    end
    
    trans_ratio_top = (v_top_dens(:,2)-cal_dens_top)./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_top.*2);
    trans_ratio_btm = (v_btm_dens(:,2)-cal_dens_btm)./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_btm.*2);
    
    trans_ratio_top_avg = trapz(phi(phi_mask),trans_ratio_top(phi_mask))./range(phi(phi_mask));
    trans_ratio_top_unc = sqrt(trapz(phi(phi_mask),(trans_ratio_top_avg-trans_ratio_top(phi_mask)).^2)./range(phi(phi_mask)));
    
    trans_ratio_btm_avg = trapz(phi(phi_mask),trans_ratio_btm(phi_mask))./range(phi(phi_mask));
    trans_ratio_btm_unc = sqrt(trapz(phi(phi_mask),(trans_ratio_btm_avg-trans_ratio_btm(phi_mask)).^2)./range(phi(phi_mask)));
    
    
    top_dens_vec(ii,:) = [top_dens_avg,top_dens_std];
    btm_ratio_vec(ii,:) = [trans_ratio_btm_avg,trans_ratio_btm_unc];
    top_ratio_vec(ii,:) = [trans_ratio_top_avg,trans_ratio_top_unc];
    out_data_vec(ii,:,:) = top_dens_norm;
    out_data_vec_ratio(ii,:) = trans_ratio_btm;
    
    
    %% Separate data into the four ports
    ports = {};
    opts_ports = {};
    %best 1. [0.274,0.95];% or 2. [0,0.95];% or [0.18,1.3];% or [0.2 1.2]
    %opts_ports.pol_lims = [0.2 1.2];%[0.274,0.95];%[0,0.95];%[0.18,1.3];%[0.18,1.3];%[0.2,1.2];%[0.15,1.1];%[0,1.5];%[0,0.95];%[0.274,0.95];%[-pi,pi];%[-0.13,pi-1.1];%[1.1,pi+0.13];%[0.8,pi];%[1.88,pi];
    [ports.top_left, ports.top_right] = separate_ports(out_data{ii}.top_halo,opts_ports);
    [ports.bottom_left, ports.bottom_right] = separate_ports(out_data{ii}.bottom_halo,opts_ports);
    
    %% Quantum correlator E
    opts_E.calc_err = do_g2_err;
    opts_E.plots = false;
    opts_E.verbose = false;
    opts_E.fit = false;
    opts_E.norm = false; %use normalised or unnormalised data
    %     opts_E.sample_proportion = 1.0;%0.1;
    lambda = 0.2;%*3;%0.25;%2.5;%10
    %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
    %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%
%     opts_E.delta_kd = [8e-3,1.5e-3,15e-3].*lambda;% main
        opts_E.delta_kd =  [4e-3,4e-3,4e-3]./4;%[0.557e-3,0.218e-3,0.782e-3];%[4e-3,4e-3,4e-3]./4;%[2e-3,1e-3,2e-3];%[4e-3,4e-3,4e-3];%[3e-3,1e-3,3e-3]./2;%[8.35e-3,1.6e-3,14e-3].*lambda; %[8.3e-3,1.6e-3,15e-3].*lambda;
    %     opts_E.delta_kd = [4.2e-3.*lambda,1.5e-3,6e-3].*4;%
    %     opts_E.delta_kd = [7.3e-3,1.2e-3,12.7e-3].*lambda;%
    opts_E.dim = 1;
    opts_E.sample_proportion = 0.001;%1.0;
    opts_E.num_samp_rep = 50;
    opts_E.do_norm = 0;
    opts_E.vol_corr = 1;
    opts_E.bin_lim = 15;
    
    %% BACK TO BACK (in the same halo)
    corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    
    global_sample_portion = 1.0;
    
    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    dkr = (dkx.*dky.*dkz).^(1/3);
    
    % variables for calculating the error
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.type='1d_vol_bb';%'1d_cart_bb';%'radial_bb';%
    corr_opts.bin_lims = 10;
    corr_opts.plots = false;
    corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
    corr_opts.fit = false;
    corr_opts.calc_err = do_g2_err;
    corr_opts.one_d_dimension = 2;
    corr_opts.two_d_dimensions = [2,3];
    %     corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*7e-3;
    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
    %     one_d_range=0.017;%0.02
    % one_d_range=0.16;
    
    one_d_range=0.015;%*3;%0.05;%0.01;%0.01;%0.017;%0.09;%0.075;%0.02%0.03
    % one_d_range=0.16;
    
    num_pts_cart = round(one_d_range./opts_E.delta_kd(corr_opts.one_d_dimension));
    num_pts_rad = round(one_d_range./dkr);
    num_pts_1 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(1)));
    num_pts_2 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(2)));
    corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
    corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);
    corr_opts.two_d_edges = {linspace(-one_d_range,one_d_range,num_pts_1)',linspace(-one_d_range,one_d_range,num_pts_2)'};
    corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';
    
    %     corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));%100 or 80 or 85 or 95
    %     corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
    
    corr_opts.rad_smoothing=nan;
    corr_opts.one_d_smoothing=nan;
    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;
    
    corr_opts.norm_samp_factor=1500;%1500;
    corr_opts.sample_proportion=0.0001;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    
    corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
    
    %% TOP HALO BACK TO BACK
    if do_g2
        
        
        [E_val(ii), corrs.ports] = E(ports,opts_E);
        if do_E_err
            E_val_ci(:,ii) = E_ci(ports,opts_E);
        end
        
        out_corrs{ii} = corrs.ports;
        
        %Expected amplitude
        %                 g12 = corrs.ports.g12.norm_g2.g2_amp(1);
        %                 g14 = corrs.ports.g14.norm_g2.g2_amp(1);
        %                 g23 = corrs.ports.g23.norm_g2.g2_amp(1);
        %                 g34 = corrs.ports.g34.norm_g2.g2_amp(1);
        corr_density='one_d_corr_density';
        if opts_E.do_norm
            if opts_E.vol_corr
                mid_pt = 1;
            else
            mid_pt=ceil(length(corrs.ports.g12.norm_g2.g2_amp)/2);
            end
            g12 = corrs.ports.g12.norm_g2.g2_amp(mid_pt);
            g14 = corrs.ports.g14.norm_g2.g2_amp(mid_pt);
            g23 = corrs.ports.g23.norm_g2.g2_amp(mid_pt);
            g34 = corrs.ports.g34.norm_g2.g2_amp(mid_pt);
        else
            if opts_E.vol_corr
                mid_pt = 1;
            else
            mid_pt=ceil(length(corrs.ports.g12.in_shot_corr.(corr_density))/2);
            end
            g12 = corrs.ports.g12.in_shot_corr.(corr_density)(mid_pt);
            g14 = corrs.ports.g14.in_shot_corr.(corr_density)(mid_pt);
            g23 = corrs.ports.g23.in_shot_corr.(corr_density)(mid_pt);
            g34 = corrs.ports.g34.in_shot_corr.(corr_density)(mid_pt);
        end
        
        
        
        top_corr_bb_vec(ii) = g12;
        btm_corr_bb_vec(ii) = g34;
        btw_1_corr_bb_vec(ii) = g23;
        btw_2_corr_bb_vec(ii) = g14;
        if opts_E.vol_corr
            g12 = corrs.ports.g12.in_shot_corr.(corr_density)(:);
            g14 = corrs.ports.g14.in_shot_corr.(corr_density)(:);
            g23 = corrs.ports.g23.in_shot_corr.(corr_density)(:);
            g34 = corrs.ports.g34.in_shot_corr.(corr_density)(:);
            E_lambda(:,ii) = -(g14+g23-g12-g34)./(g14+g23+g12+g34);
        end
        out_corrs{ii}.top_halo.corr_bb=calc_any_g2_type(corr_opts,out_data{ii}.top_halo.counts_vel');
        corr_opts.fig=['bottom halo bb corr ',num2str(phi_vec(ii))];
        out_corrs{ii}.bottom_halo.corr_bb=calc_any_g2_type(corr_opts,out_data{ii}.bottom_halo.counts_vel');
both_halo_counts = [out_data{ii}.top_halo.counts_vel';out_data{ii}.bottom_halo.counts_vel'];
corr_opts.fig=['between halo bb corr ',num2str(phi_vec(ii))];
out_corrs{ii}.between_halos.corr_bb=calc_any_g2_type(corr_opts,both_halo_counts);



top_corr_bb_vec(ii) = out_corrs{ii}.top_halo.corr_bb.in_shot_corr.one_d_corr_density(1);
        btm_corr_bb_vec(ii) = out_corrs{ii}.bottom_halo.corr_bb.in_shot_corr.one_d_corr_density(1);
        btw_1_corr_bb_vec(ii) = out_corrs{ii}.between_halos.corr_bb.in_shot_corr.one_d_corr_density(1);
        btw_2_corr_bb_vec(ii) = out_corrs{ii}.between_halos.corr_bb.in_shot_corr.one_d_corr_density(1);

        top_corr_bb_raw(:,ii) = out_corrs{ii}.top_halo.corr_bb.in_shot_corr.one_d_corr_density(:);
        btm_corr_bb_raw(:,ii) = out_corrs{ii}.bottom_halo.corr_bb.in_shot_corr.one_d_corr_density(:);
        btw_corr_bb_raw(:,ii) = out_corrs{ii}.between_halos.corr_bb.in_shot_corr.one_d_corr_density(:);

        top_corr_bb_norm(:,ii) = out_corrs{ii}.top_halo.corr_bb.norm_g2.g2_amp(:);
        btm_corr_bb_norm(:,ii) = out_corrs{ii}.bottom_halo.corr_bb.norm_g2.g2_amp(:);
        btw_corr_bb_norm(:,ii) = out_corrs{ii}.between_halos.corr_bb.norm_g2.g2_amp(:);

%                 out_corrs{ii}.top_halo.corr_bb=calc_any_g2_type(corr_opts,out_data{ii}.top_halo.counts_vel');
        %                 mid_pt=ceil(length(out_corrs{ii}.top_halo.corr_bb.norm_g2.g2_amp)/2);
        %                 top_corr_bb_vec_full(ii) = out_corrs{ii}.top_halo.corr_bb.norm_g2.g2_amp(mid_pt);
        %
        %                 corr_opts.fig=['bottom halo bb corr ',num2str(phi_vec(ii))];
        %                 out_corrs{ii}.bottom_halo.corr_bb=calc_any_g2_type(corr_opts,out_data{ii}.bottom_halo.counts_vel');
        %                 btm_corr_bb_vec_full(ii) = out_corrs{ii}.bottom_halo.corr_bb.norm_g2.g2_amp(mid_pt);
        %
        %                 both_halo_counts = [out_data{ii}.top_halo.counts_vel';out_data{ii}.bottom_halo.counts_vel'];
        %                 corr_opts.fig=['between halo bb corr',num2str(phi_vec(ii))];
        %                 out_corrs{ii}.between_halo.corr_bb=calc_any_g2_type(corr_opts,both_halo_counts);
        %                 mid_pt_bt=ceil(length(out_corrs{ii}.between_halo.corr_bb.norm_g2.g2_amp)/2);
        %                 between_corr_bb_vec_full(ii) = out_corrs{ii}.between_halo.corr_bb.norm_g2.g2_amp(mid_pt_bt);
        
        if do_g2_err
            %             top_corr_bb_unc(ii) = out_corrs{ii}.top_halo.corr_bb.norm_g2.g2_unc(1);
            %             top_corr_bb_unc(ii) = out_corrs{ii}.g14.norm_g2.g2_unc(1);
            if opts_E.vol_corr
                top_corr_bb_unc(:,ii) = out_corrs{ii}.g12.in_shot_corr.corr_unc(:);
                btm_corr_bb_vec_unc(:,ii) = out_corrs{ii}.g34.in_shot_corr.corr_unc(:);
                btw_1_corr_bb_vec_unc(:,ii) = out_corrs{ii}.g23.in_shot_corr.corr_unc(:);
                btw_2_corr_bb_vec_unc(:,ii) = out_corrs{ii}.g14.in_shot_corr.corr_unc(:);
            else
                top_corr_bb_unc(ii) = out_corrs{ii}.g12.norm_g2.g2_unc(mid_pt);
                btm_corr_bb_vec_unc(ii) = out_corrs{ii}.g34.norm_g2.g2_unc(mid_pt);
                btw_1_corr_bb_vec_unc(ii) = out_corrs{ii}.g23.norm_g2.g2_unc(mid_pt);
                btw_2_corr_bb_vec_unc(ii) = out_corrs{ii}.g14.norm_g2.g2_unc(mid_pt);
            end
        end
    end
end
%%
if do_g2
    direction_label = 'r';
    gs = {'g14','g23','g12','g34'};
    %     centers = 'rad_centers';
    %     corr_density = 'rad_corr_density';
    
    corr_density='one_d_corr_density';
    centers='x_centers';
    for ii = 1:4
        gx=gs{ii};
        stfig([gx,' comp']);
        clf
        for jj = 1:length(out_corrs)
            if ~opts_E.vol_corr
                subplot(1,3,1)
            end
            hold on
            plot(out_corrs{jj}.(gx).in_shot_corr.(centers),out_corrs{jj}.(gx).in_shot_corr.(corr_density))
            ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
            xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
            if ~opts_E.vol_corr
                subplot(1,3,2)
                hold on
                plot(out_corrs{jj}.(gx).between_shot_corr.(centers),out_corrs{jj}.(gx).between_shot_corr.(corr_density))
                ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
                xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
                subplot(1,3,3)
                hold on
                plot(out_corrs{jj}.(gx).norm_g2.(centers),out_corrs{jj}.(gx).norm_g2.g2_amp)
                ylabel(sprintf('$g^{(2)}(\\Delta %s)$',direction_label))
                xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
                
                g2_mean.(gx).val(jj)=nanmean(out_corrs{jj}.(gx).norm_g2.g2_amp);
            end
            if opts_E.vol_corr
                mid_pt = 1;
            else
                mid_pt=ceil(length(out_corrs{jj}.(gx).in_shot_corr.(corr_density))/2);
            end
            g2_raw.(gx).val(jj) = out_corrs{jj}.(gx).in_shot_corr.(corr_density)(mid_pt);
            g2_raw.(gx).vec(jj,:) = out_corrs{jj}.(gx).in_shot_corr.(corr_density)(:);
            if do_g2_err
                g2_raw.(gx).ci(jj,:,:) = out_corrs{jj}.(gx).in_shot_corr.corr_ci;
            end
            g2_raw.(gx).poisson_unc(jj,:) = 1./sqrt(cellfun(@(x) size(x,1), out_corrs{jj}.(gx).in_shot_corr.shots));
            if opts_E.do_norm
                g2_norm.(gx).vec(jj,:) = out_corrs{jj}.(gx).norm_g2.g2_amp(:);
                g2_norm.(gx).unc(jj,:) = out_corrs{jj}.(gx).norm_g2.g2_unc(:);
            end
        end
    end
end


%% Transfer percentage against phase
x = phi_vec;
y = top_dens_vec(:,1);%top_ratio_vec(:,1);%1-btm_ratio_vec(:,1);%
w = top_dens_vec(:,2);

% y = top_corr_bb_vec;%btm_ratio_vec(:,1);%
% w = top_corr_bb_unc;

% range(y)
% range(btm_ratio_vec(:,1))

% range(btm_ratio_vec(:,1))./(min(1-btm_ratio_vec(:,1))+max(1-btm_ratio_vec(:,1)))
% range(y)./(min(y)+max(y))

[x_vec, y_vec] = combine_data(phi_vec,out_data_vec(:,:,2));
% [x_vec, y_vec] = combine_data(phi_vec,1-out_data_vec_ratio);

stfig('dens over sphere');
[a, b] = sort(x_vec);
surf(a,phi,y_vec(b,:)')
% pcolor(phi,x_vec,y_vec)
% caxis([0.28 0.6])
caxis([0.25 0.9])
xlim([min(phi_vec),max(phi_vec)])
ylim([-0.2 0.2].*pi)
zlim([0 1])
ylabel('$\phi$')
xlabel('phase')
shading flat
box on
set(gca,'FontSize',19)
colorbar

[x_vec, y_vec] = combine_data(phi_vec,out_data_vec(:,:,1));
% [x_vec, y_vec] = combine_data(phi_vec,1-out_data_vec_ratio);

stfig('dens around sphere');
[a, b] = sort(x_vec);
surf(a,theta,y_vec(b,:)')
% pcolor(phi,x_vec,y_vec)
% caxis([0.28 0.6])
caxis([0.25 0.9])
xlim([min(phi_vec),max(phi_vec)])
ylim([-pi pi])
zlim([0 1])
ylabel('$\theta$')
xlabel('phase')
shading flat
box on
set(gca,'FontSize',19)
colorbar
% set(gcf,'color','w');

% xp = linspace(0,max(phi_vec));
xp = linspace(min(phi_vec),max(phi_vec));
% fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(6)).*(cos(x.*b(5) + 2*pi/b(3))) + b(4);    % Function to fit
% fit = @(b,x)  b(1).*cos(x.*b(2) + b(3)) + b(4);    % Function to fit [1.229,1,0.8088,0.906]
fit = @(b,x)  b(1).*cos(x + b(2)) + b(3);    % Function to fit
% [1.229,0.8088,0.906]
best_fit = fitnlm(x,y,fit,[0.5,0.8088,0.5],'CoefficientNames',{'Amp','Phase','Offset'}); %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
[ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'

stfig('normalised density of top halo against phase');
% clf
hold on
if plot_fit
    plot(xp,ysamp_val,'r','LineWidth',1.5)
    drawnow
    yl=ylim*1.1;
    plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
end
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
% errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
%     'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5)
% scatter(top_dens_ind(:,1),top_dens_ind(:,2))
xlabel('$\phi$')
ylabel('normalised density top halo')
grid
box on
set(gca,'FontSize',19)
xlim([min(phi_vec),max(phi_vec)])
ylim([0 1])

pred_dens = predict(best_fit,top_dens_ind(:,1));
fit_res = top_dens_ind(:,2)-pred_dens;

vis_str = string_value_with_unc(abs(2.*best_fit.Coefficients.Estimate(1)),2.*best_fit.Coefficients.SE(1),'type','b');
% vis_str = string_value_with_unc(best_fit.Coefficients.Estimate(1)./best_fit.Coefficients.Estimate(3),best_fit.Coefficients.SE(1)+best_fit.Coefficients.SE(3),'type','b');
% fprintf('\n Normalised density data\n\n')
% fprintf('\n Vis = %s\n\n',vis_str)

% %% Parameters over theta
% 
% [x_vec, y_vec] = combine_data(phi_vec,out_data_vec(:,:,1));
% 
% s_width = 10;
% m_width = 4;
% theta_c = 1:s_width;
% tot_num = floor((151-s_width)/m_width)+1;
% amp_vec = zeros(tot_num,2);
% phase_vec = zeros(tot_num,2);
% offset_vec = zeros(tot_num,2);
% theta_vec = zeros(tot_num,1);
% for ii = 1:tot_num
%     y=nanmean(y_vec(:,theta_c,1),2);
%     best_fit_c = fitnlm(x_vec,y,fit,[0.5,0.8088,0.5],'CoefficientNames',{'Amp','Phase','Offset'});
%     theta_vec(ii) = nanmean(theta(theta_c));
%     theta_c = theta_c+m_width;
%     amp_vec(ii,1) = best_fit_c.Coefficients.Estimate(1);
%     amp_vec(ii,2) = best_fit_c.Coefficients.SE(1);
%     phase_vec(ii,1) = best_fit_c.Coefficients.Estimate(2);
%     phase_vec(ii,2) = best_fit_c.Coefficients.SE(2);
%     offset_vec(ii,1) = best_fit_c.Coefficients.Estimate(3);
%     offset_vec(ii,2) = best_fit_c.Coefficients.SE(3);
% end
% %% Theta plot
% stfig('params around halo');
% % clf
% subplot(3,1,1)
% plot(theta_vec,phase_vec(:,1),'k-','LineWidth',1.5)
% hold on
% plot(theta_vec,phase_vec(:,1)-phase_vec(:,2),'r-','LineWidth',1.5)
% plot(theta_vec,phase_vec(:,1)+phase_vec(:,2),'r-','LineWidth',1.5)
% xlim([min(theta_vec),max(theta_vec)])
% xlabel('$\theta$')
% ylabel('phase')
% set(gca,'FontSize',14)
% subplot(3,1,2)
% plot(theta_vec,amp_vec(:,1),'k-','LineWidth',1.5)
% hold on
% plot(theta_vec,amp_vec(:,1)-amp_vec(:,2),'r-','LineWidth',1.5)
% plot(theta_vec,amp_vec(:,1)+amp_vec(:,2),'r-','LineWidth',1.5)
% xlim([min(theta_vec),max(theta_vec)])
% xlabel('$\theta$')
% ylabel('Amp')
% set(gca,'FontSize',14)
% subplot(3,1,3)
% plot(theta_vec,offset_vec(:,1),'k-','LineWidth',1.5)
% hold on
% plot(theta_vec,offset_vec(:,1)-offset_vec(:,2),'r-','LineWidth',1.5)
% plot(theta_vec,offset_vec(:,1)+offset_vec(:,2),'r-','LineWidth',1.5)
% xlim([min(theta_vec),max(theta_vec)])
% xlabel('$\theta$')
% ylabel('Offset')
% set(gca,'FontSize',14)
% 
% 
% %% Parameters over azimuthal angle
% 
% [x_vec, y_vec] = combine_data(phi_vec,out_data_vec(:,:,2));
% 
% s_width = 4;
% m_width = 1;
% theta_c = 1:s_width;
% tot_num = floor((151-s_width)/m_width)+1;
% amp_vec_2 = zeros(tot_num,2);
% phase_vec_2 = zeros(tot_num,2);
% offset_vec_2 = zeros(tot_num,2);
% theta_vec_2 = zeros(tot_num,1);
% for ii = 1:tot_num
%     y=nanmean(y_vec(:,theta_c),2);
%     if sum(isnan(y),1)>0
%         amp_vec_2(ii,1) = nan;
%         amp_vec_2(ii,2) = nan;
%         phase_vec_2(ii,1) = nan;
%         phase_vec_2(ii,2) = nan;
%         offset_vec_2(ii,1) = nan;
%         offset_vec_2(ii,2) = nan;
%     else
%         best_fit_c = fitnlm(x_vec,y,fit,[0.5,0.8088,0.5],'CoefficientNames',{'Amp','Phase','Offset'});
%         amp_vec_2(ii,1) = abs(best_fit_c.Coefficients.Estimate(1));
%         amp_vec_2(ii,2) = best_fit_c.Coefficients.SE(1);
%         phase_vec_2(ii,1) = best_fit_c.Coefficients.Estimate(2);
%         phase_vec_2(ii,2) = best_fit_c.Coefficients.SE(2);
%         offset_vec_2(ii,1) = best_fit_c.Coefficients.Estimate(3);
%         offset_vec_2(ii,2) = best_fit_c.Coefficients.SE(3);
%     end
%     az_vec(ii) = nanmean(phi(theta_c));
%     theta_c = theta_c+m_width;
% end
% az_mask = az_vec>-0.5 & az_vec<0.5;
% amp_vec_2_temp = amp_vec_2(az_mask,2);
% az_temp = az_vec(az_mask);
% [val,indx] = max(amp_vec_2(az_mask,1));
% max_vis_str = string_value_with_unc(val.*2,2.*amp_vec_2_temp(indx),'type','b');
% % fprintf('\n Max Vis = %s, at phi=%g\n\n',max_vis_str,az_temp(indx))
% %% Azimuthal plot
% stfig('params over halo');
% % clf
% subplot(3,1,1)
% plot(az_vec,phase_vec_2(:,1),'k-','LineWidth',1.5)
% hold on
% plot(az_vec,phase_vec_2(:,1)-phase_vec_2(:,2),'r-','LineWidth',1.5)
% plot(az_vec,phase_vec_2(:,1)+phase_vec_2(:,2),'r-','LineWidth',1.5)
% xlim([min(az_vec),max(az_vec)])
% xlabel('$\phi$')
% ylabel('phase')
% set(gca,'FontSize',14)
% xlim([-0.5,0.5])
% subplot(3,1,2)
% plot(az_vec,amp_vec_2(:,1),'k-','LineWidth',1.5)
% hold on
% plot(az_vec,amp_vec_2(:,1)-amp_vec_2(:,2),'r-','LineWidth',1.5)
% plot(az_vec,amp_vec_2(:,1)+amp_vec_2(:,2),'r-','LineWidth',1.5)
% xlim([min(az_vec),max(az_vec)])
% xlabel('$\phi$')
% ylabel('Amp')
% set(gca,'FontSize',14)
% xlim([-0.5,0.5])
% subplot(3,1,3)
% plot(az_vec,offset_vec_2(:,1),'k-','LineWidth',1.5)
% hold on
% plot(az_vec,offset_vec_2(:,1)-offset_vec_2(:,2),'r-','LineWidth',1.5)
% plot(az_vec,offset_vec_2(:,1)+offset_vec_2(:,2),'r-','LineWidth',1.5)
% xlim([min(az_vec),max(az_vec)])
% xlabel('$\phi$')
% ylabel('Offset')
% set(gca,'FontSize',14)
% xlim([-0.5,0.5])

%% g2 (and E) plots
num_shots_phi = cellfun(@(x) size(x.top_halo.num_counts,1),out_data);
ref_indx = [1];%base line points (ie points that weren't put through an interferometer)
int_indx = setdiff(1:numel(phi_vec),ref_indx);
if do_g2
    %% using the ports method
    x = phi_vec;
    y = top_corr_bb_vec;%
    l_indx = 3;
    %5

    %     E_raw=((g2_raw.g12.val(:)+g2_raw.g34.val(:)-g2_raw.g14.val(:)-g2_raw.g23.val(:))./(g2_raw.g12.val(:)+g2_raw.g34.val(:)+g2_raw.g14.val(:)+g2_raw.g23.val(:))).';
    E_raw=((g2_raw.g12.vec(:,l_indx)+g2_raw.g34.vec(:,l_indx)-1.0.*g2_raw.g14.vec(:,l_indx)-1.0.*g2_raw.g23.vec(:,l_indx))./...
        (g2_raw.g12.vec(:,l_indx)+g2_raw.g34.vec(:,l_indx)+1.0.*g2_raw.g14.vec(:,l_indx)+1.0.*g2_raw.g23.vec(:,l_indx))).';
    numerator_E = (g2_raw.g12.vec(:,l_indx)+g2_raw.g34.vec(:,l_indx)-g2_raw.g14.vec(:,l_indx)-g2_raw.g23.vec(:,l_indx)).';
    s_E = (g2_raw.g12.vec(:,l_indx)+g2_raw.g34.vec(:,l_indx)).';
    x_E = (-g2_raw.g14.vec(:,l_indx)-g2_raw.g23.vec(:,l_indx)).';
    denominator_E = (g2_raw.g12.vec(:,l_indx)+g2_raw.g34.vec(:,l_indx)+g2_raw.g14.vec(:,l_indx)+g2_raw.g23.vec(:,l_indx)).';
    %     y = btm_corr_bb_vec;%
    if do_g2_err
        if opts_E.vol_corr
            if opts_E.do_norm
                wt = g2_norm.g12.unc(:,l_indx).'./2;
            wb = g2_norm.g34.vec(:,l_indx).'./2;
            wbt1 = g2_norm.g23.vec(:,l_indx).'./2;
            wbt2 = g2_norm.g14.vec(:,l_indx).'./2;
            err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
            E_err = abs(E_raw).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
            else
            wt = top_corr_bb_unc(l_indx,:)./2;
            wb = btm_corr_bb_vec_unc(l_indx,:)./2;
            wbt1 = btw_1_corr_bb_vec_unc(l_indx,:)./2;
            wbt2 = btw_2_corr_bb_vec_unc(l_indx,:)./2;
            
            wt_n = (g2_raw.g12.vec(:,l_indx)-g2_raw.g12.ci(:,1,l_indx))./1.96;
            wt_p = (-g2_raw.g12.vec(:,l_indx)+g2_raw.g12.ci(:,2,l_indx))./1.96;
            wb_n = (g2_raw.g34.vec(:,l_indx)-g2_raw.g34.ci(:,1,l_indx))./1.96;
            wb_p = (-g2_raw.g34.vec(:,l_indx)+g2_raw.g34.ci(:,2,l_indx))./1.96;
            wbt1_n = (g2_raw.g23.vec(:,l_indx)-g2_raw.g23.ci(:,1,l_indx))./1.96;
            wbt1_p = (-g2_raw.g23.vec(:,l_indx)+g2_raw.g23.ci(:,2,l_indx))./1.96;
            wbt2_n = (g2_raw.g14.vec(:,l_indx)-g2_raw.g14.ci(:,1,l_indx))./1.96;
            wbt2_p = (-g2_raw.g14.vec(:,l_indx)+g2_raw.g14.ci(:,2,l_indx))./1.96;
            
            err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
            E_err = abs(E_raw).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
            end

        else
            wt = top_corr_bb_unc./2;
            wb = btm_corr_bb_vec_unc./2;
            wbt1 = btw_1_corr_bb_vec_unc./2;
            wbt2 = btw_2_corr_bb_vec_unc./2;
            err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
            E_err = abs(E_raw).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
        end
    else
        w=y./20;
        %take the error from poisson statistics
        wt = g2_raw.g12.vec(:,l_indx).*g2_raw.g12.poisson_unc(:,l_indx);
        wb = g2_raw.g34.vec(:,l_indx).*g2_raw.g34.poisson_unc(:,l_indx);
        wbt1 = g2_raw.g23.vec(:,l_indx).*g2_raw.g23.poisson_unc(:,l_indx);
        wbt2 = g2_raw.g14.vec(:,l_indx).*g2_raw.g14.poisson_unc(:,l_indx);
        err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2).';
        E_err = abs(E_raw).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
        E_err_deriv = sqrt(wt.^2+wb.^2).'.*(2.*abs(x_E)./denominator_E.^2)+sqrt(wbt1.^2+wbt2.^2).'.*(2.*abs(s_E)./denominator_E.^2);
    end
    
    
    
    
    xp = linspace(0,max(phi_vec));%linspace(0,pi);%;
    % fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(6)).*(cos(x.*b(5) + 2*pi/b(3))) + b(4);    % Function to fit
    fit = @(b,x)  b(1).*cos(x + b(2)) + b(3);    % Function to fit
    best_fit = fitnlm(x,y,fit,[2,0,2],'CoefficientNames',{'Amp','Phase','Offset'}); %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
    [ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    
    stfig('g2 against phase');
    clf
    hold on
    %plot(xp,ysamp_val,'r','LineWidth',1.5)
    %     drawnow
    yl=ylim*1.1;
    %plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
    colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
    %     errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    %         'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
    
    
    if opts_E.do_norm
        ht=errorbar(x+0.1,g2_norm.g12.vec(:,l_indx),wt,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hb=errorbar(x,btm_corr_bb_vec,wb,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hb=errorbar(x,g2_norm.g34.vec(:,l_indx),wb,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hbt1=errorbar(x,btw_1_corr_bb_vec,wbt1,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hbt1=errorbar(x+0.1,g2_norm.g23.vec(:,l_indx),wbt1,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hbt2=errorbar(x,btw_2_corr_bb_vec,wbt2,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hbt2=errorbar(x,g2_norm.g14.vec(:,l_indx),wbt2,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    else
    
        if ~do_g2_err
    %ht=errorbar(x,y,wt,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    ht=errorbar(x+0.1,g2_raw.g12.vec(:,l_indx),wt,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hb=errorbar(x,btm_corr_bb_vec,wb,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hb=errorbar(x,g2_raw.g34.vec(:,l_indx),wb,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hbt1=errorbar(x,btw_1_corr_bb_vec,wbt1,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hbt1=errorbar(x+0.1,g2_raw.g23.vec(:,l_indx),wbt1,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
%     hbt2=errorbar(x,btw_2_corr_bb_vec,wbt2,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hbt2=errorbar(x,g2_raw.g14.vec(:,l_indx),wbt2,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
        else
    
    ht=errorbar(x+0.1,g2_raw.g12.vec(:,l_indx),wt_n,wt_p,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);

    hb=errorbar(x,g2_raw.g34.vec(:,l_indx),wb_n,wb_p,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
    hbt1=errorbar(x+0.1,g2_raw.g23.vec(:,l_indx),wbt1_n,wbt1_p,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    
    hbt2=errorbar(x,g2_raw.g14.vec(:,l_indx),wbt2_n,wbt2_p,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
        end
        end
    %     legend([hbt1],{'g23'})
    legend([ht hb hbt1 hbt2],{'g12', 'g34','g23','g14'})
    %     scatter(x,y,'o')
    xlabel('$\Phi$')
    ylabel('$g^{(2)}_{BB}$ top halo')
    grid
    box on
    set(gca,'FontSize',19)
    
    stfig('between halos vs same halo')
    clf
    s_g2 = g2_raw.g34.vec(int_indx,l_indx)+g2_raw.g12.vec(int_indx,l_indx);
    x_g2 = g2_raw.g23.vec(int_indx,l_indx)+g2_raw.g14.vec(int_indx,l_indx);
    if ~do_g2_err
    ht=errorbar(x(int_indx),s_g2,sqrt(wt(int_indx).^2+wb(int_indx).^2),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    x_g2_unc = sqrt(wt(int_indx).^2+wb(int_indx).^2);
    s_g2_unc = sqrt(wbt1(int_indx).^2+wbt2(int_indx).^2);
    hold on
    hbt1=errorbar(x(int_indx),x_g2,sqrt(wbt1(int_indx).^2+wbt2(int_indx).^2),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
        else
    
    ht=errorbar(x(int_indx),s_g2,sqrt(wb_n.^2+wt_n.^2),sqrt(wb_p.^2+wt_p.^2),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    x_g2_unc = mean(sqrt(wb_n(int_indx).^2+wt_n(int_indx).^2),sqrt(wb_p(int_indx).^2+wt_p(int_indx).^2));
    s_g2_unc = mean(sqrt(wbt2_n(int_indx).^2+wbt1_n(int_indx).^2),sqrt(wbt1_p(int_indx).^2+wbt2_p(int_indx).^2));
    hold on
    hbt2=errorbar(x(int_indx),x_g2,sqrt(wbt2_n.^2+wbt1_n.^2),sqrt(wbt1_p.^2+wbt2_p.^2),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    end
    xlabel('Global Phase ($\Phi_L+\Phi_R$)')
    ylabel('Integrated Pair-Correlation')% 'Expected product')
    fit_2 = @(b,x)  b(1).*cos(2.*x + b(2)) + b(3);
    fit_s = fitnlm(x(int_indx),g2_raw.g34.vec(int_indx,l_indx)+g2_raw.g12.vec(int_indx,l_indx),fit_2,[2e-2,0,2e-2],'CoefficientNames',{'Amp','Phase','Offset'})
    fit_x = fitnlm(x(int_indx),g2_raw.g14.vec(int_indx,l_indx)+g2_raw.g23.vec(int_indx,l_indx),fit_2,[2e-2,0,2e-2],'CoefficientNames',{'Amp','Phase','Offset'})
    
    fit_b = fitnlm([x(int_indx),x(int_indx)+pi./2],[s_g2;x_g2],fit_2,[2e-2,0,2e-2],'CoefficientNames',{'Amp','Phase','Offset'})
    [ysamp_val,ysamp_ci]=predict(fit_b,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    plot(xp,ysamp_val,'r','LineWidth',1.5)
    
    [ysamp_val,ysamp_ci]=predict(fit_b,xp'+pi/2,'Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    plot(xp,ysamp_val,'k','LineWidth',1.5)
    ylim([0 max(max(x_g2),max(s_g2)).*1.2])
    
    %visibilities
    [vx,ix] = max(x_g2);
    [vxl,ixl] = min(x_g2);
    [vs,is] = max(s_g2);
    [vsl,isl] = min(s_g2);
    V_x = (vx-vxl)./(vx+vxl);
    dV_x = V_x.*sqrt(x_g2_unc(ix).^2+x_g2_unc(ixl).^2).*sqrt(1./(vx-vxl).^2+1./(vx+vxl).^2);
    
    V_s = (vs-vsl)./(vs+vsl);
    dV_s = V_s.*sqrt(s_g2_unc(is).^2+s_g2_unc(isl).^2).*sqrt(1./(vs-vsl).^2+1./(vs+vsl).^2);
    string_value_with_unc(V_x,dV_x,'type','b')
    string_value_with_unc(V_s,dV_s,'type','b')
    
    V_fit=max(fit_b.Coefficients.Estimate(1)./fit_b.Coefficients.Estimate(3));
    dV_fit = V_fit.*sqrt((fit_b.Coefficients.SE(1)./fit_b.Coefficients.Estimate(1)).^2+(fit_b.Coefficients.SE(3)./fit_b.Coefficients.Estimate(3)).^2);
    
    string_value_with_unc(V_fit,dV_fit,'type','b')
    
    stfig('E against phase');
    fit = @(b,x)  b(1).*cos(2.*x + b(2))+b(3);    % Function to fit
    if do_g2_err
    best_fit_E = fitnlm(x(int_indx),E_raw(int_indx),fit,[1,0.544,0],'Weight',1./E_err,'CoefficientNames',{'Amp','Phase','Offset'}) %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
    else
        best_fit_E = fitnlm(x(int_indx),E_raw(int_indx),fit,[1,0.544,0],'CoefficientNames',{'Amp','Phase','Offset'}) %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
    end
    [ysamp_val,ysamp_ci]=predict(best_fit_E,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    
    clf
    hold on
    plot(xp,ysamp_val,'r','LineWidth',1.5)
        drawnow
        yl=ylim*1.1;
    plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
    colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
    %     errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    %         'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
    hE=errorbar(x(int_indx),E_raw(int_indx),E_err(int_indx),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    plot(xp,ones(size(xp)).*nanmean(E_raw(ref_indx)),'k--','LineWidth',1.8)
    %     hb=errorbar(x,btm_corr_bb_vec,wb,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    %     hbt1=errorbar(x,btw_1_corr_bb_vec,wbt1,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    %
    %     hbt2=errorbar(x,btw_2_corr_bb_vec,wbt2,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    %     legend([hbt1],{'g23'})
    %     legend([hE],{'E($\phi$)'})
    %     scatter(x,y,'o')
    xlabel('$\Phi$')
    ylabel('E($\Phi$)')
    grid
    box on
    set(gca,'FontSize',19)
    ylim([-1 1])
    
    %probability plot
    x_P = x_g2./(s_g2+x_g2);
    s_P = s_g2./(s_g2+x_g2);
    ds_P = s_P.*sqrt((s_g2_unc.^2+x_g2_unc.^2)./(s_g2+x_g2).^2+(s_g2_unc./s_g2).^2);
    dx_P = x_P.*sqrt((s_g2_unc.^2+x_g2_unc.^2)./(s_g2+x_g2).^2+(x_g2_unc./x_g2).^2);
    fit_sp = fitnlm(x(int_indx),s_P,fit_2,[2e-2,0,2e-2],'CoefficientNames',{'Amp','Phase','Offset'});
    fit_xp = fitnlm(x(int_indx),x_P,fit_2,[2e-2,0,2e-2],'CoefficientNames',{'Amp','Phase','Offset'});
    
    stfig('Probability of ports');
    clf
    errorbar(x(int_indx),x_P,dx_P,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5)
    hold on
    errorbar(x(int_indx),s_P,ds_P,'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5)
    
    [ysamp_val,ysamp_ci]=predict(fit_sp,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    plot(xp,ysamp_val,'r','LineWidth',1.5)
    
    [ysamp_val,ysamp_ci]=predict(fit_xp,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    plot(xp,ysamp_val,'k--','LineWidth',1.5)
   
    
    xlabel('Global Phase ($\Phi_L+\Phi_R$)')
    ylabel('Probability')
    set(gca,'FontSize',17)
    legend('$P_{1,4}=P_{2,3}$','$P_{1,2}=P_{3,4}$')
    ylim([0 1])
    xlim([0.4 3.8])
    sum(E_raw(4:6).*[-1,2,1])
   %% using the entire halo
    l_indx = 5;
    num_mask = num_shots_phi>0;
   % unnormalised correlations
    top_g2_raw = top_corr_bb_raw(l_indx,num_mask);
    btm_g2_raw = btm_corr_bb_raw(l_indx,num_mask);
    btw_g2_raw = btw_corr_bb_raw(l_indx,num_mask);
   % normalised correlations
    top_g2_norm = top_corr_bb_norm(l_indx,num_mask);
    btm_g2_norm = btm_corr_bb_norm(l_indx,num_mask);
    btw_g2_norm = btw_corr_bb_norm(l_indx,num_mask);
   % E (raw and normalised)
    numerator_E = (top_g2_raw+btm_g2_raw-btw_g2_raw-btw_g2_raw);
    denominator_E = top_g2_raw+btm_g2_raw+btw_g2_raw+btw_g2_raw;

    E_calc = (top_g2_norm+btm_g2_norm-btw_g2_norm-btw_g2_norm)./(top_g2_norm+btm_g2_norm+btw_g2_norm+btw_g2_norm);
    E_raw = numerator_E./denominator_E;

% fit = @(b,x)  b(1).*cos(4.*(x.*b(4)).^b(5) + b(2))+b(3);    % Function to fit
% best_fit_E_raw = fitnlm(phi_vec(num_mask),E_raw,fit,[1,0.544,0,1,1])

    stfig('g2 raw');
    clf
    scatter(x(num_mask),top_g2_raw,50,'kx')
    hold on
    scatter(x(num_mask),btm_g2_raw,50,'rs')
    scatter(x(num_mask),btw_g2_raw,50,'bo')
    xlabel('$\Phi$')
    ylabel('$g^{(2)}_{BB}$')
    grid
    box on
    set(gca,'FontSize',19)
    legend('top', 'btm','between')

    stfig('g2 norm');
    clf
    scatter(x(num_mask),top_g2_norm,50,'kx')
    hold on
    scatter(x(num_mask),btm_g2_norm,50,'rs')
    scatter(x(num_mask),btw_g2_norm,50,'bo')
    xlabel('$\Phi$')
    ylabel('$g^{(2)}_{BB}$')
    grid
    box on
    set(gca,'FontSize',19)
    legend('top', 'btm','between')
    
    stfig('E full');
    clf
    scatter(x(num_mask),E_calc,50)
    hold on
    scatter(x(num_mask),E_raw,50)
    xlabel('$\Phi$')
    ylabel('E($\Phi$)')
    grid
    box on
    set(gca,'FontSize',19)
    ylim([-1 1])
    legend('norm', 'raw')

    %     xlim([0,max(phi_vec)])
    %     ylim([1 4])
end

%% PLOT G2'S AROUND HALOS

corr_opts.num_theta_bins = 1;
corr_opts.phi_max = 5*pi/18;
corr_opts.num_phi_bins = 5;

for ii = 1:length(phi_vec)
sph_corrs{ii}.top = spherical_g2(corr_opts,out_data{ii}.top_halo.counts_vel');
sph_corrs{ii}.btm = spherical_g2(corr_opts,out_data{ii}.bottom_halo.counts_vel');
both_halo_counts = [out_data{ii}.top_halo.counts_vel';out_data{ii}.bottom_halo.counts_vel'];
sph_corrs{ii}.btw = spherical_g2(corr_opts,both_halo_counts);
end
% 0.4 rad center for 7pi/8
%%
l_indx = 5;
phi_idx = 5;

g2_sph.top = cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),sph_corrs{phi_idx}.top.data);%,'UniformOutput',false);
g2_sph.btm = cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),sph_corrs{phi_idx}.btm.data);%,'UniformOutput',false);
g2_sph.btw = cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),sph_corrs{phi_idx}.btw.data);%,'UniformOutput',false);

stfig('g2 vs ang');
clf
plot(sph_corrs{phi_idx}.top.theta_mid,g2_sph.top)
hold on
plot(sph_corrs{phi_idx}.top.theta_mid,g2_sph.btm)
plot(sph_corrs{phi_idx}.top.theta_mid,g2_sph.btw,'k--')

stfig('E vs ang');
clf
SG=g2_sph.top+g2_sph.btm;
XG=2.*g2_sph.btw;
if corr_opts.num_phi_bins == 1
    plot(sph_corrs{phi_idx}.top.theta_mid,(SG-XG)./(SG+XG))
    xlabel('$\theta$')
elseif corr_opts.num_theta_bins == 1
    plot(sph_corrs{phi_idx}.top.phi_mid,(SG-XG)./(SG+XG))
    xlabel('$\phi$')
else
surf(sph_corrs{phi_idx}.top.phi_mid,sph_corrs{phi_idx}.top.theta_mid,(SG-XG)./(SG+XG))
colorbar
xlabel('$\phi$')
ylabel('$\theta$')
end

%%
l_indx = 5;

if corr_opts.num_phi_bins == 1
g2_sph_top = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.top.data),sph_corrs,'UniformOutput',false));
g2_sph_btm = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.btm.data),sph_corrs,'UniformOutput',false));
g2_sph_btw = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.btw.data),sph_corrs,'UniformOutput',false));

elseif corr_opts.num_theta_bins == 1
g2_sph_top = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.top.data),sph_corrs,'UniformOutput',false).');
g2_sph_btm = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.btm.data),sph_corrs,'UniformOutput',false).');
g2_sph_btw = cell2mat(cellfun(@(y) cellfun(@(x) x.in_shot_corr.one_d_corr_density(l_indx),y.btw.data),sph_corrs,'UniformOutput',false).');
end

SG=g2_sph_top+g2_sph_btm;
XG=2.*g2_sph_btw;
E_sph = ((SG-XG)./(SG+XG));
stfig('E and angle');
clf
if corr_opts.num_phi_bins == 1
surf(phi_vec(int_indx),sph_corrs{1}.top.theta_mid,E_sph(:,int_indx))
xlabel('$\Phi$')
ylabel('$\theta$')
elseif corr_opts.num_theta_bins == 1
    surf(sph_corrs{1}.top.phi_mid,phi_vec(int_indx),E_sph(int_indx,:))
    ylabel('$\Phi$ (phase)')
xlabel('$\phi$')
end

%%
corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    corr_opts.plots = false;
    corr_opts.fit = false;

corr_opts.rad_smoothing=nan;
    corr_opts.one_d_smoothing=nan;

    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;

    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    

    % variables for calculating the error
    corr_opts.calc_err = do_g2_err;
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.type='1d_vol_bb';%'1d_cart_bb';%'radial_bb';%
%%

% Halo N lims
halo_N_lims = [-1,Inf]; %2

% Num lims
data_N_lims = [3.52e3,4.66e5];%[0,Inf];

%theta bins
corr_opts.theta_bins = [0.1585,0.789];%[0,0.1]; %2
corr_opts.theta_mask_polarity = 0;

%phi bins
L_ang = 1.124;%1; % in units of 10 degrees
Ang_cen = 2.37;%1.5; % in units of 10 degrees
phi_bins_S = [-L_ang,L_ang;Ang_cen-L_ang,Ang_cen+L_ang;-L_ang,L_ang].*pi/18; %2

%bin size
opts_E.delta_kd =  [0.557e-3,0.2178e-3,0.7816e-3];%[4e-3,4e-3,4e-3]./4; %3

    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    
    
    corr_opts.bin_lims = 15;

    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];

    corr_opts.sample_proportion=0.0001;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')

S_indx = [4,5,6];

clear sph_S_corrs
for ii = 1:length(S_indx)

N_mask = out_data{S_indx(ii)}.bottom_halo.num_counts<halo_N_lims(2) & ...
out_data{S_indx(ii)}.top_halo.num_counts<halo_N_lims(2) & ...
out_data{S_indx(ii)}.bottom_halo.num_counts>halo_N_lims(1) & ...
out_data{S_indx(ii)}.top_halo.num_counts>halo_N_lims(1);

atom_num = cellfun(@(x) size(x,1), out_data{S_indx(ii)}.top_halo.counts_txy);
atom_mask = atom_num>data_N_lims(1) & atom_num<data_N_lims(2);

corr_opts.phi_bins = phi_bins_S(ii,:);

current_data.top_halo = struct_mask(out_data{S_indx(ii)}.top_halo,N_mask);
current_data.bottom_halo = struct_mask(out_data{S_indx(ii)}.bottom_halo,N_mask);

sph_S_corrs{ii}.top = spherical_section_g2(corr_opts,current_data.top_halo.counts_vel');
sph_S_corrs{ii}.btm = spherical_section_g2(corr_opts,current_data.bottom_halo.counts_vel');
both_halo_counts = [current_data.top_halo.counts_vel';current_data.bottom_halo.counts_vel'];
sph_S_corrs{ii}.btw = spherical_section_g2(corr_opts,both_halo_counts);
end

g2_S_top = cell2mat(cellfun(@(x) x.top.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));
g2_S_btm = cell2mat(cellfun(@(x) x.btm.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));
g2_S_btw = cell2mat(cellfun(@(x) x.btw.data.in_shot_corr.one_d_corr_density(:),sph_S_corrs,'UniformOutput',false));

SG=g2_S_top+g2_S_btm;
XG=2.*g2_S_btw;

N_shots_top = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.top.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));
N_shots_btm = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.btm.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));
N_shots_btw = cell2mat(cellfun(@(y) cellfun(@(x) size(x,1),y.btw.data.in_shot_corr.shots).',sph_S_corrs,'UniformOutput',false));

SG_unc = sqrt(g2_S_top./N_shots_top+g2_S_btm./N_shots_btm);
XG_unc = 2.*sqrt(g2_S_btw./N_shots_btw);

E_sph = ((SG-XG)./(SG+XG));
E_unc = E_sph.*sqrt((SG_unc.^2+XG_unc.^2)./(SG-XG).^2+(SG_unc.^2+XG_unc.^2)./(SG+XG).^2);

S_sph = sum(E_sph.*[-1,2,1],2);
S_unc = sqrt(sum((E_unc.*[1,2,1]).^2,2));

figure(11)
errorbar(1:15,S_sph,S_unc)

%% halo num
stfig('halo nums')
clf
for kk = 4;%1:length(phi_vec)
    subplot(1,2,1)
    hold on
plot(out_data{kk}.bottom_halo.num_counts)
subplot(1,2,2)
hold on
plot(out_data{kk}.top_halo.num_counts)

avg_N_btm(kk) = nanmean(out_data{kk}.bottom_halo.num_counts);
avg_N_top(kk) = nanmean(out_data{kk}.top_halo.num_counts);
end
%% Plot of E vs atom number cut off
% l_indx = 1;
% % phi_indx = 6;%1;
% N_lim_vec = 3:20;
% clear g12_vec g14_vec g23_vec g34_vec E_vec_N
% for jj = 1:length(phi_vec)
%     phi_indx = jj;
%     g12_val = out_corrs{phi_indx}.g12.in_shot_corr.shots{l_indx};
% g14_val = out_corrs{phi_indx}.g14.in_shot_corr.shots{l_indx};
% g23_val = out_corrs{phi_indx}.g23.in_shot_corr.shots{l_indx};
% g34_val = out_corrs{phi_indx}.g34.in_shot_corr.shots{l_indx};
% 
% for ii = 1:length(N_lim_vec)
%     N_lim_c = N_lim_vec(ii);
% g12_mask = out_data{phi_indx}.top_halo.num_counts(g12_val(:,1))<N_lim_c;
% g34_mask = out_data{phi_indx}.bottom_halo.num_counts(g34_val(:,1))<N_lim_c;
% 
% g14_mask = (out_data{phi_indx}.top_halo.num_counts(g14_val(:,1))+out_data{phi_indx}.bottom_halo.num_counts(g14_val(:,1)))<2.*N_lim_c;
% g23_mask = (out_data{phi_indx}.bottom_halo.num_counts(g23_val(:,1))+out_data{phi_indx}.top_halo.num_counts(g23_val(:,1)))<2.*N_lim_c;
% 
% g12_vec(jj,ii) = sum(out_corrs{phi_indx}.g12.in_shot_corr.shots{l_indx}(g12_mask,3));
% g14_vec(jj,ii) = sum(out_corrs{phi_indx}.g14.in_shot_corr.shots{l_indx}(g14_mask,3));
% g23_vec(jj,ii) = sum(out_corrs{phi_indx}.g23.in_shot_corr.shots{l_indx}(g23_mask,3));
% g34_vec(jj,ii) = sum(out_corrs{phi_indx}.g34.in_shot_corr.shots{l_indx}(g34_mask,3));
% 
% E_vec_N(jj,ii) = (g12_vec(jj,ii)+g34_vec(jj,ii)-g23_vec(jj,ii)-g14_vec(jj,ii))./(g12_vec(jj,ii)+g34_vec(jj,ii)+g23_vec(jj,ii)+g14_vec(jj,ii));
% end
% end
%%
function settings_list = settings_list_func();
settings_list = {
        'c:\remote\settings202122Jul092216.xml';%evap: ~0.838 MHz
        'c:\remote\settings202122Jul092456.xml';%evap: ~0.8385 MHz
        'c:\remote\settings202120Jul151102.xml';%evap: ~0.839 MHz
        'c:\remote\settings202122Jul092533.xml';%evap: ~0.8395 MHz
        'c:\remote\settings202120Jul151039.xml';%evap: ~0.840 MHz
        'c:\remote\settings202122Jul092648.xml';%evap: ~0.8405 MHz
        'c:\remote\settings202120Jul151012.xml';%evap: ~0.841 MHz
        'c:\remote\settings202120Jul150821.xml';%evap: ~0.8415 MHz
        'c:\remote\settings202120Jul150601.xml';%evap: ~0.842 MHz
        'c:\remote\settings202120Jul150519.xml';%evap: ~0.8425 MHz
        'c:\remote\settings202120Jul150504.xml';%evap: ~0.843 MHz
        'c:\remote\settings202120Jul150351.xml';%evap: ~0.8435 MHz
        'c:\remote\settings202120Jul150327.xml';%evap: ~0.844 MHz
        'c:\remote\settings202120Jul150310.xml';%evap: ~0.8445 MHz
        'c:\remote\settings202120Jul142856.xml';%evap: ~0.845 MHz
        'c:\remote\settings202120Jul151131.xml';%evap: ~0.846 MHz
        'c:\remote\settings202120Jul151208.xml';%evap: ~0.847 MHz
        'c:\remote\settings202122Jul161245.xml';%evap: ~0.848 MHz
        'c:\remote\settings202122Jul161222.xml';%evap: ~0.849 MHz
        'c:\remote\settings202122Jul161154.xml';%evap: ~0.850 MHz
        'c:\remote\settings202122Jul161425.xml';%evap: ~0.851 MHz
        'c:\remote\settings202122Jul161039.xml';%evap: ~0.852 MHz
        'c:\remote\settings202122Jul161117.xml';%evap: ~0.853 MHz
        'c:\remote\settings202122Jul161444.xml';%evap: ~0.854 MHz
        };
end
%%
function out=cell_comp(a,b)
if size(a,1)~=size(b,1) 
    out= false; 
else
    out= all(cellfun(@isequal,a,b));
end
end