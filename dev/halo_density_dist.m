%% Initializing path
clear all;
% close all;
this_folder = fileparts(which(mfilename));%'C:\Users\helium\cloudstor\PROJECTS\';
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
%% Import directory
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\pulse_characterisation\';%'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';
% data_folder='20210727_k=0,-1_halo_mirror_tests\98_kHz_2';%1.227641917170454e-01
% data_folder='20210728_k=0,-1_halo_splitter_tests\test_1';
% data_folder='20210728_k=0,-1_halo_splitter_tests\test_2';
% data_folder='98_kHz_3';%2.874773680074205e-01
% data_folder='98_kHz_4';%1.868659019918052e-01
% data_folder='98_kHz_5';%2.384441983254298e-01
% data_folder='98_kHz';%1.914509594759020e-01
% data_folder='102_kHz';%1.242533737685670e-01
% data_folder='104_kHz';%1.450573430876497e-01
% data_folder='112_kHz';%3.387780300875264e-01
% data_folder='wider_pulse_1';%2.188036094952446e-01
% data_folder='wider_pulse_2';%1.589455442811602e-01
% data_folder='wider_pulse_3';%1.663914007678803e-01
% data_folder='wider_pulse_4';%1.857936635571151e-01
% data_folder='wider_pulse_5';%2.733980619224182e-01
% data_folder='20210728_k=0,-1_halo_splitter_tests\test_6';
% data_folder='wider_pulse_6';%1.345175790431930e-01
% data_folder='k=0,-1,-2_halos_data\tighter_trap\20210325_k=0,-1,-2_halos_8';
% data_folder='k=0,-1,-2_halos_data\20210302_k=0,-1,-2_halo_data_4'; 
% data_folder='20210808_k=0,-1,-2_rt_scan_mid_trap_equal_delay_6'; 
% data_folder='20210803_k=0,-1,-2_rt_scan_mid_trap_3';
% data_folder={'20210727_k=0,-1_halo_mirror_tests\98_kHz_2','20210727_k=0,-1_halo_mirror_tests\98_kHz','20210727_k=0,-1_halo_mirror_tests\98_kHz_3',...
%     '20210727_k=0,-1_halo_mirror_tests\98_kHz_4','20210727_k=0,-1_halo_mirror_tests\98_kHz_5'};
%%'20210728_k=0,-1_halo_splitter_tests\test_4';%test_4
%%'20210727_k=0,-1_halo_mirror_tests\98_kHz_2';%
% 'pulse_characterisation\20210201_bragg_pulse_testing\prototyping\gauss_2';
% 'pulse_characterisation\20210129_mirror_test_1\sqrt(20)';
data_folder = '20230404_k=0,-1_and_mirror_Vsh_0_75_Vq_0_7';%'test data for halo_density_dist'; %'test data for halo_density_dist';%'20230405_check_splitter_50_50';%20230322_k=+1,0_and_splitter_2_500_mus
% data_folder='20210715_testing_bs_on_k=0,-1_halo\moving_detuning_further';
% data_folder='20210715_testing_bs_on_k=0,-1_halo\moving_detuning';
% data_folder='20210715_testing_bs_on_k=0,-1_halo\first_pulse_change_detuning';
% data_folder='20210715_testing_bs_on_k=0,-1_halo\first_pulse_5';
% data_folder = '20210712_k=0,-1,-2_halos_new_trap';
% data_folder = '20210713_k=0,-1,-2_halos_new_trap_2';
% data_folder='20210325_k=0,-1,-2_halos_8';
% data_folder='20210325_k=0,-1_mirror_test_tight_trap\gaussian_test';
%data_folder='20210219_mirror_800_mus_chirped';
% data_folder='20210330_mz_scan_tighter_trap_5_mus_amp_12.9';
%data_folder='20210325_k=0,-1_mirror_test_tight_trap\13_111_kHz';
% data_folder='k=0,-1,-2_halos_data\20210301_k=0,-1,-2_halo_data_3';
% data_folder='k=0,-1,-2_halos_data\20210302_k=0,-1,-2_halo_data_4';
% data_folder='k=0,-1,-2_halos_data\tighter_trap\20210325_k=0,-1,-2_halos_8';
% data_folder='20210223_k=0,-1_norm_evap_848';
% data_folder='20210217_k=0,-1_norm_evap_755';
% data_folder='20210218_bragg_pulses_at_1200_mus\splitter_3';
% data_folder='20210312_bragg_scan\1200_mus_delay\8_111_kHz';
% data_folder='20210310_bragg_pulse_test\800mus_mirror_v5';
% data_folder='20230301_second_attempt_k=0,-1_mj=0_only_halo_Vq_2_Vsh_0.65';
% data_folder='20200901_k=0,-1_transfer_vs_amp\Pamp_16';%'20210201_bragg_pulse_testing\scan\amp_5_5';
% data_folder='20211206_scaning_across_freq\norm_0_844_MHz';
% data_folder='20211027_validating_bragg_pulses\bs_2';
opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;

% opts.import.shot_num = 4:28; %

%validating first run
% opts.import.shot_num = 32:41; %0 can select specific shots to import
% opts.import.shot_num = 42:51; %0.035 can select specific shots to import
% opts.import.shot_num = 52:58; %0.045 can select specific shots to import
% opts.import.shot_num = 62:71; %0.03 can select specific shots to import
% opts.import.shot_num = 72:81; %0.02 can select specific shots to import
% opts.import.shot_num = 82:101; %0.01 can select specific shots to import
% opts.import.shot_num = 103:196; %0.025 can select specific shots to import

%validating second run
% opts.import.shot_num = 16:26; %0.025 and 1:16 also
% opts.import.shot_num = 27:40; %0.015 and 1:16 also
% opts.import.shot_num = 41:52; %0.018 and 1:16 also
% opts.import.shot_num = 53:61; %0.018 and 1:16 also

%validate bs
% opts.import.shot_num = 171:179;
% opts.import.shot_num = 180:188;

% opts.import.shot_num = 199:249;%0
% opts.import.shot_num = 250:262;%+0.02
% opts.import.shot_num = 263:275;%+0.01
% opts.import.shot_num = 276:286;%+0.008
% opts.import.shot_num = 287:327;%+0.006

%% Chose which halo(s) to analyse
opts.do_top_halo = 1;% analyse the top halo?
opts.do_btm_halo = 1;% analyse the bottom halo?

%% Chose if you want to look at a narrow or wide slice of the halo
slice_type = 'medium';
if strcmp(slice_type,'narrow')
    opts.vel_conv.top.z_mask = [-0.4,0.4];%
    opts.vel_conv.btm.z_mask = [-0.4,0.4];%in units of radius ([-0.68,0.68])
elseif strcmp(slice_type,'extra wide')
    opts.vel_conv.top.z_mask = [-0.87,0.87];%
    opts.vel_conv.btm.z_mask = [-0.87,0.87];%in units of radius ([-0.68,0.68])
elseif strcmp(slice_type,'medium')
    opts.vel_conv.top.z_mask = [-0.6,0.6];%
    opts.vel_conv.btm.z_mask = [-0.6,0.6];%in units of radius ([-0.68,0.68])
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
% 
opts.vel_conv.top.z_mask = [-0.9,0.9];%[-1,1];%[-0.9,0.9];
opts.vel_conv.btm.z_mask = [-0.9,0.9];%[-1,1];%[-0.9,0.9];%in units of radius ([-0.68,0.68])

radius_lim = [0.056,0.072];%[0.00,0.14];%[0.03,0.09];%[0.01.*0.065,0.08];%[0.01.*0.065,0.1];%[0.05,0.07];%[0.3,1.61].*0.065;%[0.61,1.26].*0.065;%[0.9,1.05];%[0.89,1.11];
ang_lim = 90;%angular limit in degrees

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,6];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 0e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = 0;%2;%10;%0;% %minimum allowed number in halo 10

opts.plot_dist = false; %do you want to see all the detailed stuff about the halo distributions

opts.cent.nan_cull = false; %do you want to cull nan's

%% Background stuff
if iscell(data_folder)
    cli_header('Setting up for %s and', data_folder{1});
    for ii = 2:length(data_folder)
        cli_header(' %s ', data_folder{ii});
    end
else
    cli_header('Setting up for %s', data_folder);
end
opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
opts.data_src = fullfile(opts.data_root, data_folder);
opts.data_dir = data_folder;
opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
% opts.index.filename = sprintf('index__%s__%.0f', opts.data_dir);
opts.label = data_folder;
opts.tag = 0;
opts.full_out = false;
opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds
opts.shot_bounds = [];
combined_struct = @(S,T) cell2struct(cellfun(@vert_or_horz_cat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);
% if ~exist(opts.fig_dir, 'dir')
%     mkdir(opts.fig_dir);
% end
% Run the function!

%% Set up out dir
%set up an output dir %https://gist.github.com/ferryzhou/2269380
% if (exist([opts.data_src, '\out'], 'dir') == 0), mkdir(fullfile(opts.data_src, '\out')); end
% %make a subfolder with the ISO timestamp for that date
% anal_out.dir = sprintf('%sout\\%s\\', ...
%     [opts.data_src, '\'], datestr(datetime('now'), 'yyyymmddTHHMMSS'));
% if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end

%% import raw data
dirs_list = fullfile(opts.data_root, data_folder);
if iscell(dirs_list)
    for ii = 1:length(dirs_list)
        opts.import.dir = dirs_list{ii};
        opts.import.cache_save_dir = fullfile(dirs_list{ii}, 'cache', 'import\');
        
        if ii == 1
            [data, ~] = import_mcp_tdc_data(opts.import);
        else
            [data_temp, ~] = import_mcp_tdc_data(opts.import);
            data = combined_struct(data,data_temp);
        end
    end
else
    [data, ~] = import_mcp_tdc_data(opts.import);
end
% [data, ~] = import_mcp_tdc_data(opts.import);

%% remove any ringing
% data_ht_spot=hotspot_mask(data);
% data.counts_txy=data_ht_spot.masked.counts_txy;
% data.num_counts=data_ht_spot.masked.num_counts;
opts.ring_lim = 0.09e-6;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
data_masked = ring_removal(data,opts.ring_lim);

%% set up relevant constants
hebec_constants

%% find centers
opts.cent.visual = 0; %from 0 to 2
opts.cent.savefigs = 0;
opts.cent.correction = 0;
opts.cent.correction_opts.plots = 0;

opts.cent.top.visual = 0; %from 0 to 2
opts.cent.top.savefigs = 0;
opts.cent.top.threshold = [130,5000,5000].*1e3;%130
opts.cent.top.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[16,7,10].*1e3;
opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.top.method = {'margin','average','average'};

opts.cent.mid.visual = 0; %from 0 to 2
opts.cent.mid.savefigs = 0;
opts.cent.mid.threshold = [130,5000,5000].*1e3;
opts.cent.mid.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[16,7,10].*1e3;
opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.mid.method = {'margin','average','average'};

opts.cent.btm.visual = 0; %from 0 to 2
opts.cent.btm.savefigs = 0;
opts.cent.btm.threshold = [130,5000,5000].*1e3;%[130,2000,2000].*1e3;
opts.cent.btm.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.btm.method = {'margin','average','average'};

% opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% opts.cent.t_bounds = {[1.741,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% opts.cent.t_bounds = {[2.134,2.148],[2.148,2.161],[2.161,2.18],[2.13,2.2]};
 opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.884,3.896],[3.75,4]};%
%  opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
% opts.cent.t_bounds = {[5.350,5.356],[5.361,5.367],[5.372,5.380],[5.34,5.39]};%time bounds for the different momentum states (for full evap settings)
bec = halo_cent(data_masked,opts.cent);

%% run some checks
% atoms number
% laser maybe?
num_check = data_masked.num_counts>opts.num_lim;
% num_masked = data_masked.num_counts;
% num_masked(~num_check) = NaN;
% num_outlier = isoutlier(num_masked);
% ~num_outlier &
is_shot_good = num_check & bec.centre_OK_mid';
if opts.do_top_halo
    is_shot_good = is_shot_good & bec.centre_OK_top';
end
if opts.do_btm_halo
    is_shot_good = is_shot_good & bec.centre_OK_btm';
end
data_masked_halo = struct_mask(data_masked,is_shot_good);
bec_masked_halo = struct_mask(bec,is_shot_good);

%% Find the velocity widths
opts.bec_width.g0 = const.g0;
opts.bec_width.fall_time = 0.417;
bec_masked_halo = bec_width_txy_to_vel(bec_masked_halo,opts.bec_width);

%% convert data to velocity
% zero velocity point
t0 = bec_masked_halo.centre_top(:,1);%ones(size(bec_masked_halo.centre_top,1),1).*1.76;%;%%.*2.1749;%.*3.8772;%72;%
x0 = bec_masked_halo.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
y0 = bec_masked_halo.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;

%% generate top halo
opts.vel_conv.top.visual = 0;
opts.vel_conv.top.plot_percentage = 0.95;
opts.vel_conv.top.title = 'top halo';
opts.vel_conv.top.const.g0 = const.g0;
opts.vel_conv.top.const.fall_distance = const.fall_distance;
opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.top.v_mask=radius_lim; %bounds on radisu in units of m/s
opts.vel_conv.top.ang_lim = ang_lim; %angular limits of the azimuthal angle
opts.vel_conv.top.y_mask = [-1.9,1.9]; %in units of radius
opts.vel_conv.top.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point

opts.vel_conv.top.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s

opts.vel_conv.top.bec_center.north = bec_masked_halo.centre_top;
opts.vel_conv.top.bec_center.south = bec_masked_halo.centre_mid;
opts.vel_conv.top.bec_width.north = bec_masked_halo.width_top;
opts.vel_conv.top.bec_width.south = bec_masked_halo.width_mid;

%%
if opts.do_top_halo
    top_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.top);
else
    top_halo_intial.counts_vel = {};
    top_halo_intial.counts_vel_norm = {};
    top_halo_intial.num_counts = [];
end

%% generate bottom halo
opts.vel_conv.btm.visual = 0;
opts.vel_conv.btm.plot_percentage = 0.5;
opts.vel_conv.btm.title = 'bottom halo';
opts.vel_conv.btm.const.g0 = const.g0;
opts.vel_conv.btm.const.fall_distance = const.fall_distance;
opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.btm.v_mask=radius_lim;%[0.61,1.26];%[0.3,1.61];%[0.61,1.26];%[0.89,1.11]; %bounds on radisu as multiple of radius value
opts.vel_conv.btm.ang_lim = ang_lim; %angular limits of the azimuthal angle
opts.vel_conv.btm.y_mask = [-1.9,1.9]; %in units of radius
opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point

opts.vel_conv.btm.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s

opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;

%%
if opts.do_btm_halo
    bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);
else
    bottom_halo_intial.counts_vel = {};
    bottom_halo_intial.counts_vel_norm = {};
    bottom_halo_intial.num_counts = [];
end

%% mask out halos with nums to low
halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim;
halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim;
if opts.do_top_halo && opts.do_btm_halo
    halo_N_check = halo_N_check_top & halo_N_check_btm;
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

%% plot some histogram checks

%% Setting up variables
halos.top_halo = top_halo;
halos.bottom_halo = bottom_halo;
halos.bec = bec_halo;
opts.plot_opts = [];
opts.plot_opts.only_dists = true;
opts.plot_opts.const.g0 = const.g0;
opts.plot_opts.const.fall_distance = const.fall_distance;
% if opts.plot_dist
%     plot_checks(halos,opts.plot_opts);
% end

d = opts.plot_opts.const.fall_distance;
g0 = opts.plot_opts.const.g0;
tf = sqrt(2*d/g0);%arrival time of zero velocity particles

v_top_zxy = cell2mat(top_halo.counts_vel_norm);
v_top_zxy_unnorm = cell2mat(top_halo.counts_vel);
r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
r_dist_top_unnorm = sqrt(v_top_zxy_unnorm(:,1).^2+v_top_zxy_unnorm(:,2).^2+v_top_zxy_unnorm(:,3).^2);
N_top = top_halo.num_counts;

v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm);
v_btm_zxy_unnorm = cell2mat(bottom_halo.counts_vel);
if ~isempty(v_btm_zxy)
    r_dist_btm = sqrt(v_btm_zxy(:,1).^2+v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2);
    r_dist_btm_unnorm = sqrt(v_btm_zxy_unnorm(:,1).^2+v_btm_zxy_unnorm(:,2).^2+v_btm_zxy_unnorm(:,3).^2);
else
    r_dist_btm = [];
    r_dist_btm_unnorm = [];
end
N_btm = bottom_halo.num_counts;

num_shots = size(top_halo.counts_vel_norm,1);
%     top_halo = halos.top_halo;
%     bottom_halo = halos.bottom_halo;
%     bec_masked = halos.bec;
%%
% for ii = 1:size(top_halo.counts_vel,1)
%     top_halo.counts_vel
%     bottom_halo.counts_vel
% end

%% histograming
nbins=151;%51;%50;%
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
for ii = 1:(nbins-1)
    r_btm_zxy_masked = r_dist_btm_unnorm(theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1));
    r_top_zxy_masked = r_dist_top_unnorm(theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1));
    v_btm_r(ii,1) = mean(r_btm_zxy_masked);
    v_top_r(ii,1) = mean(r_top_zxy_masked);
%     v_btm_dens(ii,1) = size(r_btm_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
%     v_top_dens(ii,1) = size(r_top_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));

    r_btm_zxy_masked = r_dist_btm_unnorm(phi_bins(ii)<phi_btm & phi_btm<=phi_bins(ii+1));
    r_top_zxy_masked = r_dist_top_unnorm(phi_bins(ii)<phi_top & phi_top<=phi_bins(ii+1));
    v_btm_r(ii,2) = mean(r_btm_zxy_masked);
    v_top_r(ii,2) = mean(r_top_zxy_masked);
%     v_btm_dens(ii,2) = size(r_btm_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));
%     v_top_dens(ii,2) = size(r_top_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));

%     theta(ii) = mean(theta_bins(ii:(ii+1)));
%     phi(ii) = mean(phi_bins(ii:(ii+1)));
end
% v_btm_dens = v_btm_dens./size(bottom_halo.counts_vel,1);
% v_top_dens = v_top_dens./size(top_halo.counts_vel,1);
v_btm_dens = [];
v_top_dens = [];
v_btm_dens_unc = [];
v_top_dens_unc = [];
phi_mask_top = (phi_top<0.154& phi_top>-0.154);
phi_mask_btm = (phi_btm<0.154& phi_btm>-0.154);
phi_sig = 0.05;
r_btm_zxy_masked=smooth_hist(theta_btm(phi_mask_btm),'sigma',0.1,'lims',[-pi,pi],'bin_num',nbins);
r_top_zxy_masked=smooth_hist(theta_top(phi_mask_top),'sigma',0.1,'lims',[-pi,pi],'bin_num',nbins);
v_btm_dens(:,1) = r_btm_zxy_masked.count_rate.smooth./num_shots;
v_top_dens(:,1) = r_top_zxy_masked.count_rate.smooth./num_shots;
v_btm_dens_unc(:,1) = sqrt(r_btm_zxy_masked.count_rate.smooth).*sqrt(abs(r_btm_zxy_masked.bin.edge(1:end-1)...
    -r_btm_zxy_masked.bin.edge(2:end)));
v_top_dens_unc(:,1) = sqrt(r_top_zxy_masked.count_rate.smooth).*sqrt(abs(r_top_zxy_masked.bin.edge(1:end-1)...
    -r_top_zxy_masked.bin.edge(2:end)));


r_btm_zxy_masked=smooth_hist(phi_btm,'sigma',phi_sig,'lims',[-pi/2,pi/2],'bin_num',nbins);
r_top_zxy_masked=smooth_hist(phi_top,'sigma',phi_sig,'lims',[-pi/2,pi/2],'bin_num',nbins);
v_btm_dens(:,2) = r_btm_zxy_masked.count_rate.smooth./num_shots;
v_top_dens(:,2) = r_top_zxy_masked.count_rate.smooth./num_shots;

v_top_dens_2d = hist3([theta_top phi_top],'Nbins',[nbins nbins]);
v_btm_dens_2d = hist3([theta_btm phi_btm],'Nbins',[nbins nbins]);

theta = linspace(-pi,pi,nbins);
phi = linspace(-pi/2,pi/2,nbins);

%% 2D projections in velocity space

t=linspace(0,2*pi,1e4);
if opts.do_top_halo
    stfig('density of top halo');
    clf
    subplot(2,3,1)
    ndhist(v_top_zxy(:,2:3));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,2)
    ndhist(v_top_zxy(:,[1,3]));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,3)
    ndhist(v_top_zxy(:,1:2));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
    
    subplot(2,3,4)
    ndhist(v_top_zxy_unnorm(:,2:3));
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,5)
    ndhist(v_top_zxy_unnorm(:,[1,3]));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,6)
    ndhist(v_top_zxy_unnorm(:,1:2));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
end

if opts.do_btm_halo
    stfig('density of bottom halo');
    clf
    subplot(2,3,1)
    ndhist(v_btm_zxy(:,2:3));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,2)
    ndhist(v_btm_zxy(:,[1,3]));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,3)
    ndhist(v_btm_zxy(:,1:2));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
    
    subplot(2,3,4)
    ndhist(v_btm_zxy_unnorm(:,2:3));
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,5)
    ndhist(v_btm_zxy_unnorm(:,[1,3]));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,6)
    ndhist(v_btm_zxy_unnorm(:,1:2));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
end


%% Average radius in spherical coordinates
stfig('average radius vs angle');
clf
subplot(2,1,1)
plot(theta(2:end),v_btm_r(:,1),'linewidth',1.5)
hold on
plot(theta(2:end),v_top_r(:,1),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\theta\)')
subplot(2,1,2)
plot(phi(2:end),v_btm_r(:,2),'linewidth',1.5)
hold on
plot(phi(2:end),v_top_r(:,2),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\phi\)')

%% Density in spherical coordinates
stfig('density vs angle');
clf
subplot(2,1,1)
plot(theta./pi,v_btm_dens(:,1),'linewidth',1.5)
hold on
plot(theta./pi,v_top_dens(:,1),'linewidth',1.5)
xlim([-1,1])
legend('bottom','top')
ylabel('Average density')
xlabel('\(\theta/\pi\)')
subplot(2,1,2)
plot(phi./pi,v_btm_dens(:,2),'linewidth',1.5)
hold on
plot(phi./pi,v_top_dens(:,2),'linewidth',1.5)
xlim([asin(opts.vel_conv.top.z_mask(1)),asin(opts.vel_conv.top.z_mask(2))]./pi)
% ylim([0 6000])
legend('bottom','top')
ylabel('Average density')
xlabel('\(\phi/\pi\)')

%% Density against height in z
stfig('density vs height');
% rad_mean = mean(bottom_halo_intial.rad);
rad_mean = 1;
clf
plot(rad_mean.*sin(phi),v_btm_dens(:,2),'linewidth',1.5)
hold on
plot(rad_mean.*sin(phi),v_top_dens(:,2),'linewidth',1.5)
legend('bottom','top')
ylabel('Average density')
xlabel('\(v_z\) (m/s)')

%% ratio in spherical coordinates
if opts.do_btm_halo && opts.do_top_halo
    stfig('density ratio vs angle');
%     clf
    subplot(2,1,1)
    theta_ratio = (v_btm_dens(:,1))./(v_top_dens(:,1));
    plot(theta./pi,theta_ratio,'linewidth',1.5)
    v_ratio_unc = theta_ratio.*sqrt((v_btm_dens_unc./v_btm_dens(:,1)).^2+(v_top_dens_unc./v_top_dens(:,1)).^2);
    hold on
%     plot(theta./pi,(v_btm_dens(:,1))./(v_top_dens(:,1))+v_ratio_unc,...
%         theta./pi,(v_btm_dens(:,1))./(v_top_dens(:,1))-v_ratio_unc,'r-','linewidth',1.5)
    xlim([-1,1])
    ylabel('density ratio')
    xlabel('\(\theta/\pi\)')
    subplot(2,1,2)
    hold on
    plot(phi./pi,(v_btm_dens(:,2))./(v_top_dens(:,2)),'linewidth',1.5)
    ylabel('density ratio')
    xlabel('\(\phi/\pi\)')
    
    xlim([asin(opts.vel_conv.top.z_mask(1)),asin(opts.vel_conv.top.z_mask(2))]./pi)
%     xlim([-0.1,0.1])
%     if strcmp(slice_type,'narrow')
%         xlim([-0.153,0.153])
%     else
%         xlim([-1,1]./2)
%     end
    
%     stfig('desnity ratio 2d');
%     clf
%     pcolor(theta./pi,phi./pi,v_btm_dens_2d./v_top_dens_2d)
%     xlabel('\(\theta\)')
%     ylabel('\(\phi\)')
%     shading flat
%     colorbar
%     caxis([0.1 1.1])
    top_dens_norm = (v_top_dens(:,2))./(v_btm_dens(:,2)+v_top_dens(:,2));
%     temp_norm = (temp_top(:,2))./(temp_btm(:,2)+temp_top(:,2));
    
%     trans_ratio_top = (top_dens_norm-1+temp_norm)./(1-2.*(1-temp_norm));
%     trans_ratio_top = (v_top_dens(:,2)-cal_dens_top(:,2))./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_top(:,2).*2);
%     trans_ratio_btm = (v_btm_dens(:,2)-cal_dens_btm(:,2))./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_btm(:,2).*2);
    
    stfig('density visibility vs angle');
%     clf
    grid on
    box on
    
hold on
    plot(phi./pi,1-top_dens_norm,'linewidth',1.5)
    ylabel('density vis')
    xlabel('\(\phi/\pi\)')
    xlim([-1,1])
%     xlim([asin(opts.vel_conv.top.z_mask(1)),asin(opts.vel_conv.top.z_mask(2))]./pi)
    phi_mask = phi<0.1 & phi>-0.1;
    top_dens_avg = trapz(phi(phi_mask),top_dens_norm(phi_mask))./range(phi(phi_mask))
    top_dens_std = sqrt(trapz(phi(phi_mask),(top_dens_avg-top_dens_norm(phi_mask)).^2)./range(phi(phi_mask)))
end

%% density plot in full spherical coordinates
stfig('spherical density plot');
for ii = 1:(nbins)
    for jj = 1:(nbins)
        ang_mask_btm = theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1) ...
            & phi_bins(jj)<phi_btm & phi_btm<=phi_bins(jj+1);
        area = abs(theta_bins(ii+1)-theta_bins(ii))*abs(sin(phi_bins(jj+1))-sin(phi_bins(jj)));
        ang_mask_top = theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1) ...
            & phi_bins(jj)<phi_top & phi_top<=phi_bins(jj+1);
        density_btm(ii,jj) = sum(ang_mask_btm)./area;
        avg_r_btm(ii,jj)  = nanmean(r_dist_btm_unnorm(ang_mask_btm));
        density_top(ii,jj) = sum(ang_mask_top)./area;
        avg_r_top(ii,jj)  = nanmean(r_dist_top_unnorm(ang_mask_top));
    end
end
subplot(2,2,1)
pcolor(phi,theta,avg_r_btm)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('average radius bottom')
subplot(2,2,2)
title('average density bottom')
pcolor(phi,theta,density_btm)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('density')
subplot(2,2,3)
pcolor(phi,theta,avg_r_top)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('average radius top')
subplot(2,2,4)
title('average density top')
pcolor(phi,theta,density_top)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('density')

%% 2d comparison
%mean(top_halo.rad)
%mean(bottom_halo.rad)
rad_shift=0.059;
z_shift_top = [rad_shift.*ones(size(v_top_zxy_unnorm,1),1),zeros(size(v_top_zxy_unnorm,1),2)];
z_shift_btm = [rad_shift.*ones(size(v_btm_zxy_unnorm,1),1),zeros(size(v_btm_zxy_unnorm,1),2)];
if opts.do_top_halo && opts.do_btm_halo
    stfig('density of halos');
    clf
    subplot(1,2,1)
    combined_vzxy = [v_top_zxy_unnorm+z_shift_top;...
        v_btm_zxy_unnorm-z_shift_btm];
    ndhist(combined_vzxy(:,[3,1]),'bins',5);
    xlabel('$v_y$ (m/s)')
    ylabel('$v_z$ (m/s)')
    subplot(1,2,2)
    ndhist(combined_vzxy(:,[2,1]),'bins',5);
    xlabel('$v_x$ (m/s)')
    ylabel('$v_z$ (m/s)')
    colormap('default')
end

%% 2d comparison radius
if opts.do_top_halo && opts.do_btm_halo
    stfig('density of halos radius');
    clf
    combined_vzr = [v_top_zxy_unnorm(:,1)+z_shift_top(:,1),rxy_top;...
        v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),rxy_btm];
    ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
    xlabel('$v_r$ (m/s)')
    ylabel('$v_z$ (m/s)')
    colormap('default')
    axis equal
    caxis([0 0.5e7])
end

%% 2d comparison radius
% if opts.do_top_halo && opts.do_btm_halo
%     stfig('density of halos radius');
% %     clf
%     combined_vzr = [v_top_zxy_unnorm(:,1)+z_shift_top(:,1),rxy_top;...
%         v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),rxy_btm];
%     ndhist(combined_vzr(:,[1,2]),'bins',4,'filter');
%     xlabel('$v_z$ (m/s)')
%     ylabel('$v_r$ (m/s)')
%     colormap('default')
%     axis equal
%     caxis([0 2000000])
%     xlim([-0.115 0.115])
%     ylim([0.005 0.08])
%     set(gca,'FontSize',18)
%     
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% end
%% full 3d comparison
plt_p = 1;
stfig('halo comparison');
clf
plot_mask_top = rand(size(v_top_zxy_unnorm,1),1)<plt_p;
plot_mask_btm = rand(size(v_btm_zxy_unnorm,1),1)<plt_p;
scatter3(v_top_zxy_unnorm(plot_mask_top,2),v_top_zxy_unnorm(plot_mask_top,3),v_top_zxy_unnorm(plot_mask_top,1)+z_shift_top(plot_mask_top,1),'r.')
hold on
scatter3(v_btm_zxy_unnorm(plot_mask_btm,2),v_btm_zxy_unnorm(plot_mask_btm,3),v_btm_zxy_unnorm(plot_mask_btm,1)-z_shift_btm(plot_mask_btm,1),'b.')
xlabel('$v_z$')
ylabel('$v_y$')
zlabel('$v_z$')
axis equal
%%
stfig('radial distribution');
% clf
r_hist_top=smooth_hist(r_dist_top,'sigma',0.0001);
r_hist_btm=smooth_hist(r_dist_btm,'sigma',0.0001);
r_hist_top_un=smooth_hist(r_dist_top_unnorm,'sigma',0.0001);
r_hist_btm_un=smooth_hist(r_dist_btm_unnorm,'sigma',0.0001);
subplot(2,1,1)
plot(r_hist_top.bin.centers,r_hist_top.counts.smooth,'linewidth',1.5)
hold on
plot(r_hist_btm.bin.centers,r_hist_btm.counts.smooth,'linewidth',1.5)
xlabel('r')
ylabel('Freq')
xlim([min([r_hist_top.bin.centers;r_hist_btm.bin.centers]),...
    max([r_hist_top.bin.centers;r_hist_btm.bin.centers])])
legend('top','btm')
subplot(2,1,2)
plot(r_hist_top_un.bin.centers,r_hist_top_un.counts.smooth,'linewidth',1.5)
hold on
plot(r_hist_btm_un.bin.centers,r_hist_btm_un.counts.smooth,'linewidth',1.5)
ylimit = max([r_hist_top_un.counts.smooth;r_hist_btm_un.counts.smooth]);
plot([0.130159/2 0.130159/2],[-0.1,ylimit.*2],'k-','linewidth',1.5)
ylim([0 ylimit.*1.1])
xlabel('r')
ylabel('Freq')
xlim([min([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers]),...
    max([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers])])
legend('top','btm','expected radius')
%     if plot_opts.only_dists
%         error('breaking out')
%     end
%
%%
stfig('distribution of halo radius');
% clf
vr_hist_top=smooth_hist(top_halo.rad,'sigma',0.00003);
vr_hist_btm=smooth_hist(bottom_halo.rad,'sigma',0.00003);
plot(vr_hist_top.bin.centers,vr_hist_top.counts.smooth,'linewidth',1.5)
hold on
plot(vr_hist_btm.bin.centers,vr_hist_btm.counts.smooth,'linewidth',1.5)
xlabel('radius of halo')
ylabel('Freq')
legend('top','btm')