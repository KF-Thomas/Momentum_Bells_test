%% E vs angle

%% Initializing path
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
    %
%     '20210803_k=0,-1,-2_rt_scan_mid_trap_equal_delay_1'
%     '20210804_k=0,-1,-2_rt_scan_mid_trap_equal_delay_2'%g34 not biased downward, good E
%     %     '20210805_k=0,-1,-2_rt_scan_mid_trap_equal_delay_3'% g34 phi=pi baised downward very strongly (not good E)
%     '20210806_k=0,-1,-2_rt_scan_mid_trap_equal_delay_4'% g34 phi=pi baised downward
%     '20210807_k=0,-1,-2_rt_scan_mid_trap_equal_delay_5'% g34 phi=pi baised downward but only slightly
%     '20210808_k=0,-1,-2_rt_scan_mid_trap_equal_delay_6'% g34 phi=pi baised downward
    
%     '20211028_k=0,-1,-2_RT_run_1'
%     '20211028_k=0,-1,-2_RT_run_2'
%     '20211029_k=0,-1,-2_RT_run_3'
%     '20211029_k=0,-1,-2_RT_run_4'
%     
%     '20211030_k=0,-1,-2_RT_run_5'
%     '20211031_k=0,-1,-2_RT_run_6'
    
    '20211117_adj_k=0,-1,-2_RT_run_23' %[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,3*pi/4,7*pi/8,pi]+1.053./2
    '20211118_adj_k=0,-1,-2_RT_run_24'
    '20211118_adj_k=0,-1,-2_RT_run_25'
    '20211119_adj_k=0,-1,-2_RT_run_26'
    '20211119_adj_k=0,-1,-2_RT_run_27'
    '20211120_adj_k=0,-1,-2_RT_run_28'
    '20211122_adj_k=0,-1,-2_RT_run_29'
    '20211123_adj_k=0,-1,-2_RT_run_30'
    '20211124_adj_k=0,-1,-2_RT_run_31'

    
    };

force_reimport = false;
force_reimport_norm = true;

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.5e3;%9e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper = Inf;%Inf;%20;%3.5;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_norm = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper_norm = Inf;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_both = -1;

y_cut = 11e-3;

% opts.import.shot_num = 1:16; %can select which shots you want to import

%% Calibration settings
plot_fit = true;
do_g2=true;
do_g2_err=true;
do_E_err=false;
do_range_cut=false;
out_data = {};
phi_vec = [];

num_shots_norm = 0;
opts.tag = 1;

%% Run over each folder
for folder_indx = 1:length(data_folders)
    %% import raw data
    data_folder = data_folders{folder_indx};
    opts.import.dir = fullfile(opts.data_root, data_folder);
    opts.logfile = fullfile(opts.import.dir,log_folder);
    opts.loglabfile = fullfile(opts.import.dir,log_lab_folder);
    phi_logs = table2array(readtable(opts.logfile));
    opts.cent.nan_cull = true;
    opts.import.force_reimport = force_reimport;
    opts.import.force_cache_load = ~opts.import.force_reimport;
    %         opts.import.shot_num = 1:24; %can select which shots you want to import
    
    
    %% Chose which halo(s) to analyse
    opts.do_top_halo = 1;% analyse the top halo?
    opts.do_btm_halo = 1;% analyse the bottom halo?
    
    %% Chose if you want to look at a narrow or wide slice of the halo
    opts.vel_conv.top.z_mask = [-0.9,0.9];
    opts.vel_conv.btm.z_mask = [-0.9,0.9];%in units of radius ([-0.68,0.68])
    radius_lim = [0.0,0.13];%[0.03,0.09];%[0.0,0.12];%[0.06,0.069];%[0.03,0.09];%[0.0585,0.069];%[0.06,0.067];%[0.06,0.07];%[0.045,0.085];%[0.79,1.17].*0.065;%[0.61,1.26];%[0.89,1.11];%[0.89,1.16];%[0.9,1.05];%
    ang_lim = 90;
    
    %% Import parameters
    tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
    tmp_ylim=[-35e-3, 35e-3];
    tlim=[0,4];
    opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];
    
    opts.num_lim = 2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
    
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
    
    %% import raw data
    [data, ~] = import_mcp_tdc_data(opts.import);
    
    %% remove any ringing or hotspots
    data_ht_spot=hotspot_mask(data);
    data.counts_txy=data_ht_spot.masked.counts_txy;
    data.num_counts=data_ht_spot.masked.num_counts;
    opts.ring_lim = -1;%0.09e-6;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
    data_masked = ring_removal(data,opts.ring_lim);
    
    %% add labview import
    if opts.tag
        opts_tab = detectImportOptions(opts.loglabfile);
        opts_tab.Delimiter = {','};
        logs = readtable(opts.loglabfile,opts_tab);
        tags = logs{:,6};
        %% select a specific shot type if you wish
        shot_type = 'main';
        tag_mask = cellfun(@(x) strcmp(x, shot_type), tags');
        tag_mask = [tag_mask,zeros(1,length(data_masked.num_counts)-length(tags))];
        data_masked = struct_mask(data_masked,logical(tag_mask),1);
    end
    
    
    %% set up relevant constants
    hebec_constants
    
    %% find centers
    opts.cent.visual = 0; %from 0 to 2
    opts.cent.savefigs = 0;
    opts.cent.correction = 0;
    opts.cent.correction_opts.plots = 0;
    
    opts.cent.top.visual = 0; %from 0 to 2
    opts.cent.top.savefigs = 0;
    opts.cent.top.threshold = [130,6000,6000].*1e3;
    opts.cent.top.min_threshold = [0,8,8].*1e3;%[0,3,3].*1e3;%[16,7,10].*1e3;
    opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.top.method = {'margin','average','average'};
    
    opts.cent.mid.visual = 0; %from 0 to 2
    opts.cent.mid.savefigs = 0;
    opts.cent.mid.threshold = [130,6000,6000].*1e3;
    opts.cent.mid.min_threshold = [0,8,8].*1e3;%[16,7,10].*1e3;
    opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.mid.method = {'margin','average','average'};
    
    opts.cent.btm.visual = 0; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [130,6000,6000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [0,8,8].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.btm.method = {'margin','average','average'};
    
    %              opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
    opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
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
    is_shot_good = is_shot_good & phi_log_check';
    if do_range_cut
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
    phi_logs_masked = phi_log_matched(is_shot_good,:);
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
    
    opts.vel_conv.top.centering_correction = [0,0,0];
    
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
    opts.vel_conv.btm.plot_percentage = 0.95;
    opts.vel_conv.btm.title = 'bottom halo';
    opts.vel_conv.btm.const.g0 = const.g0;
    opts.vel_conv.btm.const.fall_distance = const.fall_distance;
    opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
    opts.vel_conv.btm.v_mask=radius_lim; %bounds on radisu as multiple of radius value
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
    halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim & top_halo_intial.num_counts<opts.halo_N_lim_upper;
    halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim & bottom_halo_intial.num_counts<opts.halo_N_lim_upper;
    halo_N_check_both = (top_halo_intial.num_counts+bottom_halo_intial.num_counts)>opts.halo_N_lim_both;
    
    if opts.do_top_halo && opts.do_btm_halo
        halo_N_check = halo_N_check_top & halo_N_check_btm &halo_N_check_both;
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
    
    %% Seperate halos by phase
    
    unique_phi = unique(phi_logs_masked(halo_N_check,1)); %unique phases used in this scan
    for ii = 1:length(unique_phi)
        current_phi = unique_phi(ii,1);
        phi_mask = (phi_logs_masked(halo_N_check,1)==current_phi);
        masked_top = struct_mask(top_halo,phi_mask');
        masked_bottom = struct_mask(bottom_halo,phi_mask');
        if ~ismember(current_phi,phi_vec)
            phi_vec = [phi_vec,current_phi];
        end
        phi_indx = find(phi_vec==current_phi);
        if length(out_data)<phi_indx
            out_data{phi_indx} = {};
            out_data{phi_indx}.top_halo = masked_top;
            out_data{phi_indx}.bottom_halo = masked_bottom;
        else
            out_data{phi_indx}.top_halo = combined_struct(out_data{phi_indx}.top_halo,masked_top);
            out_data{phi_indx}.bottom_halo =combined_struct(out_data{phi_indx}.bottom_halo,masked_bottom);
        end
    end
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

%% interate over angular slice
ang_lim_vec = 5:5:25;
halo_num_lim = 10;%17;
halo_N_lim_both = 1;
rad_lim = [0.03,0.09];%

for ii = 1:length(phi_vec)
    for kk = 1:length(ang_lim_vec)
        %% Separate data into the four ports
        ang_lim = ang_lim_vec(kk);
        temp_data.top_halo.counts_vel = out_data{ii}.top_halo.counts_vel;
        temp_data.bottom_halo.counts_vel = out_data{ii}.bottom_halo.counts_vel;
        for jj = 1:size(temp_data.top_halo.counts_vel,1)
            c_counts_t = temp_data.top_halo.counts_vel{jj};
            c_counts_b = temp_data.bottom_halo.counts_vel{jj};
            ang_mask_top = abs(180/pi*atan(c_counts_t(:,1)./sqrt(c_counts_t(:,2).^2+c_counts_t(:,3).^2)))<ang_lim;
            ang_mask_btm = abs(180/pi*atan(c_counts_b(:,1)./sqrt(c_counts_b(:,2).^2+c_counts_b(:,3).^2)))<ang_lim;
            radius_mask_top = (c_counts_t(:,1).^2+c_counts_t(:,2).^2+c_counts_t(:,3).^2)<(rad_lim(2)).^2 ...
                & (c_counts_t(:,1).^2+c_counts_t(:,2).^2+c_counts_t(:,3).^2)>(rad_lim(1)).^2;
            radius_mask_btm = (c_counts_b(:,1).^2+c_counts_b(:,2).^2+c_counts_b(:,3).^2)<(rad_lim(2)).^2 ...
                & (c_counts_b(:,1).^2+c_counts_b(:,2).^2+c_counts_b(:,3).^2)>(rad_lim(1)).^2;
            mask_top = ang_mask_top & radius_mask_top;
            mask_btm = ang_mask_btm & radius_mask_btm;
            if sum(mask_top)>halo_num_lim || sum(mask_btm)>halo_num_lim || (sum(mask_top)+sum(mask_btm))<halo_N_lim_both
                mask_top = logical(zeros(size(ang_mask_top)));
                mask_btm = logical(zeros(size(ang_mask_btm)));
            end
            
            temp_data.top_halo.counts_vel{jj} = c_counts_t(mask_top,:);
            temp_data.bottom_halo.counts_vel{jj} = c_counts_b(mask_btm,:);
        end
        
        ports = {};
        opts_ports.norm = false;
        [ports.top_left, ports.top_right] = separate_ports(temp_data.top_halo,opts_ports);
        [ports.bottom_left, ports.bottom_right] = separate_ports(temp_data.bottom_halo,opts_ports);
        
        %% Quantum correlator E
        opts_E.calc_err = do_g2_err;
        opts_E.plots = true;
        opts_E.verbose = false;
        opts_E.fit = true;
        opts_E.norm = false; %use normalised or unnormalised data
        %     opts_E.sample_proportion = 1.0;%0.1;
        lambda = 0.2;%0.25;%2.5;%10
        %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
        %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%
        opts_E.delta_kd = [8e-3,1.5e-3,15e-3].*lambda;%
        %     opts_E.delta_kd = [8.4e-3,6e-3,24e-3]./4;
        %     opts_E.delta_kd = [4.2e-3.*lambda,1.5e-3,6e-3].*4;%
        %     opts_E.delta_kd = [7.3e-3,1.2e-3,12.7e-3].*lambda;%
        opts_E.dim = 1;
        opts_E.sample_proportion = 0.01;%1.0;
        opts_E.num_samp_rep = 50;
        opts_E.do_norm = 0;
        opts_E.vol_corr = 1;
        opts_E.bin_lim = 60;
        
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
        corr_opts.type='1d_cart_bb';%'radial_bb';%
        corr_opts.plots = true;
        corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
        corr_opts.fit = false;
        corr_opts.calc_err = do_g2_err;
        corr_opts.one_d_dimension = 2;
        corr_opts.two_d_dimensions = [2,3];
        %     corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*7e-3;
        corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
        %     one_d_range=0.017;%0.02
        % one_d_range=0.16;
        
        one_d_range=0.05;%0.01;%0.01;%0.017;%0.09;%0.075;%0.02%0.03
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
        corr_opts.direction_labels = {'z','x','y'};
        corr_opts.low_mem=true;
        
        corr_opts.norm_samp_factor=1500;%1500;
        corr_opts.sample_proportion=1.0;%1.0;%0.65;%1500;
        corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
        corr_opts.do_pre_mask=false;
        corr_opts.sorted_dir=nan;
        corr_opts.sort_norm=0;
        
        corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
        
        %% TOP HALO BACK TO BACK
        if do_g2
            
            
            [E_val(kk,ii), corrs.ports] = E(ports,opts_E);
            if do_E_err
                E_val_ci(:,kk,ii) = E_ci(ports,opts_E);
            end
            
            out_corrs{kk,ii} = corrs.ports;
            
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
            
            
            
            top_corr_bb_vec(kk,ii) = g12;
            btm_corr_bb_vec(kk,ii) = g34;
            btw_1_corr_bb_vec(kk,ii) = g23;
            btw_2_corr_bb_vec(kk,ii) = g14;
            if opts_E.vol_corr
                g12 = corrs.ports.g12.in_shot_corr.(corr_density)(:);
                g14 = corrs.ports.g14.in_shot_corr.(corr_density)(:);
                g23 = corrs.ports.g23.in_shot_corr.(corr_density)(:);
                g34 = corrs.ports.g34.in_shot_corr.(corr_density)(:);
                E_lambda(:,kk,ii) = -(g14+g23-g12-g34)./(g14+g23+g12+g34);
            end
            if do_g2_err
                if opts_E.vol_corr
                    top_corr_bb_unc(:,kk,ii) = out_corrs{kk,ii}.g12.in_shot_corr.corr_unc(:);
                btm_corr_bb_vec_unc(:,kk,ii) = out_corrs{kk,ii}.g34.in_shot_corr.corr_unc(:);
                btw_1_corr_bb_vec_unc(:,kk,ii) = out_corrs{kk,ii}.g23.in_shot_corr.corr_unc(:);
                btw_2_corr_bb_vec_unc(:,kk,ii) = out_corrs{kk,ii}.g14.in_shot_corr.corr_unc(:);
                else
                top_corr_bb_unc(kk,ii) = out_corrs{kk,ii}.g12.norm_g2.g2_unc(mid_pt);
                btm_corr_bb_vec_unc(kk,ii) = out_corrs{kk,ii}.g34.norm_g2.g2_unc(mid_pt);
                btw_1_corr_bb_vec_unc(kk,ii) = out_corrs{kk,ii}.g23.norm_g2.g2_unc(mid_pt);
                btw_2_corr_bb_vec_unc(kk,ii) = out_corrs{kk,ii}.g14.norm_g2.g2_unc(mid_pt);
                end
            end
        end
    end
end
%% interate over radial mask
rad_lim_vec = [];%linspace(0.01,0.05,5);
ang_lim = 25;
for ii = 1:length(phi_vec)
    for kk = 1:length(rad_lim_vec)
        %% Separate data into the four ports
        dr = rad_lim_vec(kk);
        rad_lim = 0.063+[-dr,dr];
        temp_data.top_halo.counts_vel = out_data{ii}.top_halo.counts_vel;
        temp_data.bottom_halo.counts_vel = out_data{ii}.bottom_halo.counts_vel;
        for jj = 1:size(temp_data.top_halo.counts_vel,1)
            c_counts_t = temp_data.top_halo.counts_vel{jj};
            c_counts_b = temp_data.bottom_halo.counts_vel{jj};
            ang_mask_top = abs(180/pi*atan(c_counts_t(:,1)./sqrt(c_counts_t(:,2).^2+c_counts_t(:,3).^2)))<ang_lim;
            ang_mask_btm = abs(180/pi*atan(c_counts_b(:,1)./sqrt(c_counts_b(:,2).^2+c_counts_b(:,3).^2)))<ang_lim;
            radius_mask_top = (c_counts_t(:,1).^2+c_counts_t(:,2).^2+c_counts_t(:,3).^2)<(rad_lim(2)).^2 ...
                & (c_counts_t(:,1).^2+c_counts_t(:,2).^2+c_counts_t(:,3).^2)>(rad_lim(1)).^2;
            radius_mask_btm = (c_counts_b(:,1).^2+c_counts_b(:,2).^2+c_counts_b(:,3).^2)<(rad_lim(2)).^2 ...
                & (c_counts_b(:,1).^2+c_counts_b(:,2).^2+c_counts_b(:,3).^2)>(rad_lim(1)).^2;
            mask_top = ang_mask_top & radius_mask_top;
            mask_btm = ang_mask_btm & radius_mask_btm;
            if sum(mask_top)>halo_num_lim || sum(mask_btm)>halo_num_lim
                mask_top = logical(zeros(size(ang_mask_top)));
                mask_btm = logical(zeros(size(ang_mask_btm)));
            end
            
            temp_data.top_halo.counts_vel{jj} = c_counts_t(mask_top,:);
            temp_data.bottom_halo.counts_vel{jj} = c_counts_b(mask_btm,:);
        end
        
        ports = {};
        opts_ports.norm = false;
        [ports.top_left, ports.top_right] = separate_ports(temp_data.top_halo,opts_ports);
        [ports.bottom_left, ports.bottom_right] = separate_ports(temp_data.bottom_halo,opts_ports);
        
        %% Quantum correlator E
        opts_E.calc_err = do_g2_err;
        opts_E.plots = true;
        opts_E.verbose = false;
        opts_E.fit = true;
        opts_E.norm = false; %use normalised or unnormalised data
        %     opts_E.sample_proportion = 1.0;%0.1;
        lambda = 0.2;%0.25;%2.5;%10
        %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
        %     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%
        opts_E.delta_kd = [8e-3,1.5e-3,15e-3].*lambda;%
        %     opts_E.delta_kd = [8.4e-3,6e-3,24e-3]./4;
        %     opts_E.delta_kd = [4.2e-3.*lambda,1.5e-3,6e-3].*4;%
        %     opts_E.delta_kd = [7.3e-3,1.2e-3,12.7e-3].*lambda;%
        opts_E.dim = 1;
        opts_E.sample_proportion = 0.01;%1.0;
        opts_E.num_samp_rep = 50;
        opts_E.do_norm = 0;
        opts_E.vol_corr = 1;
        opts_E.bin_lim = 60;
        
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
        corr_opts.type='1d_cart_bb';%'radial_bb';%
        corr_opts.plots = true;
        corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
        corr_opts.fit = false;
        corr_opts.calc_err = do_g2_err;
        corr_opts.one_d_dimension = 2;
        corr_opts.two_d_dimensions = [2,3];
        %     corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*7e-3;
        corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
        %     one_d_range=0.017;%0.02
        % one_d_range=0.16;
        
        one_d_range=0.05;%0.01;%0.01;%0.017;%0.09;%0.075;%0.02%0.03
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
        corr_opts.direction_labels = {'z','x','y'};
        corr_opts.low_mem=true;
        
        corr_opts.norm_samp_factor=1500;%1500;
        corr_opts.sample_proportion=1.0;%1.0;%0.65;%1500;
        corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
        corr_opts.do_pre_mask=false;
        corr_opts.sorted_dir=nan;
        corr_opts.sort_norm=0;
        
        corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
        
        %% TOP HALO BACK TO BACK
        if do_g2
            
            
            [E_val(kk,ii), corrs.ports] = E(ports,opts_E);
            if do_E_err
                E_val_ci(:,kk,ii) = E_ci(ports,opts_E);
            end
            
            out_corrs_rad{kk,ii} = corrs.ports;
            
            %Expected amplitude
            %                 g12 = corrs.ports.g12.norm_g2.g2_amp(1);
            %                 g14 = corrs.ports.g14.norm_g2.g2_amp(1);
            %                 g23 = corrs.ports.g23.norm_g2.g2_amp(1);
            %                 g34 = corrs.ports.g34.norm_g2.g2_amp(1);
            corr_density='one_d_corr_density';
            if opts_E.do_norm
                mid_pt=ceil(length(corrs.ports.g12.norm_g2.g2_amp)/2);
                g12 = corrs.ports.g12.norm_g2.g2_amp(mid_pt);
                g14 = corrs.ports.g14.norm_g2.g2_amp(mid_pt);
                g23 = corrs.ports.g23.norm_g2.g2_amp(mid_pt);
                g34 = corrs.ports.g34.norm_g2.g2_amp(mid_pt);
            else
                mid_pt=ceil(length(corrs.ports.g12.in_shot_corr.(corr_density))/2);
                g12 = corrs.ports.g12.in_shot_corr.(corr_density)(mid_pt);
                g14 = corrs.ports.g14.in_shot_corr.(corr_density)(mid_pt);
                g23 = corrs.ports.g23.in_shot_corr.(corr_density)(mid_pt);
                g34 = corrs.ports.g34.in_shot_corr.(corr_density)(mid_pt);
            end
            
            
            
            top_corr_bb_vec_rad(kk,ii) = g12;
            btm_corr_bb_vec_rad(kk,ii) = g34;
            btw_1_corr_bb_vec_rad(kk,ii) = g23;
            btw_2_corr_bb_vec_rad(kk,ii) = g14;
            if opts_E.vol_corr
                g12 = corrs.ports.g12.in_shot_corr.(corr_density)(:);
                g14 = corrs.ports.g14.in_shot_corr.(corr_density)(:);
                g23 = corrs.ports.g23.in_shot_corr.(corr_density)(:);
                g34 = corrs.ports.g34.in_shot_corr.(corr_density)(:);
                E_lambda_rad(:,kk,ii) = -(g14+g23-g12-g34)./(g14+g23+g12+g34);
            end
            if do_g2_err
                top_corr_bb_unc_rad(kk,ii) = out_corrs_rad{kk,ii}.g12.norm_g2.g2_unc(mid_pt);
                btm_corr_bb_vec_unc_rad(kk,ii) = out_corrs_rad{kk,ii}.g34.norm_g2.g2_unc(mid_pt);
                btw_1_corr_bb_vec_unc_rad(kk,ii) = out_corrs_rad{kk,ii}.g23.norm_g2.g2_unc(mid_pt);
                btw_2_corr_bb_vec_unc_rad(kk,ii) = out_corrs_rad{kk,ii}.g14.norm_g2.g2_unc(mid_pt);
            end
        end
    end
end
%%
for ii = 1:length(phi_vec)
phi_c = phi_vec(ii);
stfig(['E(',num2str(phi_c),') vs lambda and theta']);
clf
surf(ang_lim_vec,(1:opts_E.bin_lim).*lambda,E_lambda(:,:,ii))
xlabel('$\theta$ (degrees)')
ylabel('$\lambda$')
end
%%
if do_g2
    direction_label = 'r';
    gs = {'g14','g23','g12','g34'};
    %     centers = 'rad_centers';
    %     corr_density = 'rad_corr_density';
    corr_density='one_d_corr_density';
    centers='x_centers';
    wt = top_corr_bb_unc./2;
wb = btm_corr_bb_vec_unc./2;
wbt1 = btw_1_corr_bb_vec_unc./2;
wbt2 = btw_2_corr_bb_vec_unc./2;
err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
denom = zeros(size(err_tot));
    for ii = 1:4
        gx=gs{ii};
        %         stfig([gx,' comp']);
        %         clf
        for ww = 1:size(out_corrs,2)
            for zz = 1:size(out_corrs,1)
                %             if ~opts_E.vol_corr
                %                 subplot(1,3,1)
                %             end
                %             hold on
                %             plot(out_corrs{jj}.(gx).in_shot_corr.(centers),out_corrs{jj}.(gx).in_shot_corr.(corr_density))
                %             ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
                %             xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
                %             if ~opts_E.vol_corr
                %                 subplot(1,3,2)
                %                 hold on
                %                 plot(out_corrs{jj}.(gx).between_shot_corr.(centers),out_corrs{jj}.(gx).between_shot_corr.(corr_density))
                %                 ylabel(sprintf('$G^{(2)}(\\Delta %s)$ coincedence density',direction_label))
                %                 xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
                %                 subplot(1,3,3)
                %                 hold on
                %                 plot(out_corrs{jj}.(gx).norm_g2.(centers),out_corrs{jj}.(gx).norm_g2.g2_amp)
                %                 ylabel(sprintf('$g^{(2)}(\\Delta %s)$',direction_label))
                %                 xlabel(sprintf('$\\Delta %s$ Seperation',direction_label))
                %
                %                 g2_mean.(gx).val(jj)=nanmean(out_corrs{jj}.(gx).norm_g2.g2_amp);
                %             end
                if opts_E.vol_corr
                    mid_pt = 1;
                else
                    mid_pt=ceil(length(out_corrs{zz,ww}.(gx).in_shot_corr.(corr_density))/2);
                end
                try
                    g2_raw.(gx).val(zz,ww) = out_corrs{zz,ww}.(gx).in_shot_corr.(corr_density)(mid_pt);
                    denom(:,zz,ww) = denom(:,zz,ww) + out_corrs{zz,ww}.(gx).in_shot_corr.(corr_density)(:);
                catch
                    g2_raw.(gx).val(zz,ww) = nan;
                    denom(:,zz,ww) = denom(:,zz,ww) + nan;
                end
            end
        end
    end
    E_unc = sqrt(err_tot.^2.*(1+E_lambda.^2)./denom.^2);
end
%% E and g vs lambda

l_lim = 8;
l_min = 0.2;
ang_indx = 3;%the angle limit
phi_indx = 1;%which phase should we use
lambda_vec = (1:opts_E.bin_lim).*lambda;
l_mask = lambda_vec>l_min;
lambda_vec = lambda_vec(l_mask);


mdl = @(b,x) rob_E(6.42912831449739e-01.*[x,x,x]',b(1),0)' + b(3);
% mdl = @(b,x) rob_E(b(2).*[x,x,x]',b(1),0)' + b(3);
strt_pt = [10,1,0];
% w = 1./E_unc(:,ii).^2;
mdl_fit=fitnlm(lambda_vec,E_lambda(l_mask,ang_indx,phi_indx),mdl,strt_pt,'Weight',1./E_unc(l_mask,ang_indx,phi_indx).^0);
xp=linspace(0,l_lim,1000);
[ysamp_val,ysamp_ci]=predict(mdl_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));

mdl_fit_2=fitnlm(lambda_vec,E_lambda(l_mask,ang_indx,2),mdl,strt_pt,'Weight',1./E_unc(l_mask,ang_indx,2).^0.0);
[ysamp_val_2,ysamp_ci_2]=predict(mdl_fit_2,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));


stfig('E vs lambda');
clf
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7];
colors_main=colors_main./255;
if do_g2_err
errorbar(6.42912831449739e-01.*lambda_vec,E_lambda(l_mask,ang_indx,phi_indx),E_unc(l_mask,ang_indx,phi_indx),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
hold on
errorbar(6.42912831449739e-01.*lambda_vec,E_lambda(l_mask,ang_indx,2),E_unc(l_mask,ang_indx,2),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(3,:),'LineWidth',2.5)
else
scatter(lambda_vec,E_lambda(l_mask,ang_indx,phi_indx),'o',...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
end
% errorbar(lambda,E_val(:,1),E_unc(:,1),'kx')
hold on
plot(6.42912831449739e-01.*xp,ysamp_val,'r','LineWidth',1.5)
plot(6.42912831449739e-01.*xp,ysamp_val_2,'b','LineWidth',1.5)
legend('$\Phi=0$','$\Phi=\pi/2$')

%     drawnow
%     yl=ylim*1.1;
% plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
xlim([0,6.42912831449739e-01.*l_lim])
ylim([-1 1])
xlabel('$\lambda$')
ylabel(['E($\Phi$)'])
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)

stfig('g2 vs lambda');
clf
gx = 'g12';
V2=(1:60)'.^6.*prod([8e-3,1.5e-3,15e-3].*lambda).^2;
if do_g2_err
    y = out_corrs{ang_indx,phi_indx}.(gx).in_shot_corr.(corr_density)(l_mask)./V2(l_mask);
    y2 = out_corrs{ang_indx,2}.(gx).in_shot_corr.(corr_density)(l_mask)./V2(l_mask);
    errorbar(lambda_vec,y,wt(l_mask,ang_indx,phi_indx)./V2(l_mask),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
hold on
errorbar(lambda_vec,y2,wt(l_mask,ang_indx,2)./V2(l_mask),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(3,:),'LineWidth',2.5)
else
scatter(lambda_vec,E_lambda(l_mask,ang_indx,phi_indx),'o',...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
end
mdl_2 = @(b,x) rob_G2([x,x,x]',[b],0)';
mdl_fit_g2=fitnlm(lambda_vec,y2,mdl_2,[22,1,1e10]);%,'Weight',w_vec);
[ysamp_val,ysamp_ci]=predict(mdl_fit_g2,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));

xlabel('$\lambda$')
ylabel('$G^{(2)}_{12}$')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)


%% E and g vs slice
l_lim = 8;
l_min = 0.05;
l_indx =3;%the angle limit
phi_indx = 1;%which phase should we use
lambda_vec = (1:opts_E.bin_lim).*lambda;
l_mask = lambda_vec>l_min;
lambda_vec = lambda_vec(l_mask);


mdl = @(b,x) rob_E(b(2).*[x,x,x]',b(1),0)' + b(3);
% w = 1./E_unc(:,ii).^2;
% mdl_fit=fitnlm(lambda_vec,E_lambda(l_mask,ang_indx,phi_indx),mdl,[10,1,0],'Weight',1./E_unc(l_mask,ang_indx,phi_indx).^0.0);
% xp=linspace(0,l_lim,1000);
% [ysamp_val,ysamp_ci]=predict(mdl_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));
% 
% mdl_fit_2=fitnlm(lambda_vec,E_lambda(l_mask,ang_indx,2),mdl,[10,1,0],'Weight',1./E_unc(l_mask,ang_indx,2).^0.0);
% [ysamp_val_2,ysamp_ci_2]=predict(mdl_fit_2,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));


stfig('E vs ang');
clf
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7];
colors_main=colors_main./255;
if do_g2_err
errorbar(ang_lim_vec,E_lambda(l_indx,:,phi_indx),E_unc(l_indx,:,phi_indx),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
hold on
errorbar(ang_lim_vec,E_lambda(l_indx,:,2),E_unc(l_indx,:,2),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_main(3,:),'LineWidth',2.5)
else
scatter(ang_lim_vec,E_lambda(l_indx,:,phi_indx),'o',...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
end
% errorbar(lambda,E_val(:,1),E_unc(:,1),'kx')
% hold on
% plot(xp,ysamp_val,'r','LineWidth',1.5)
% plot(xp,ysamp_val_2,'b','LineWidth',1.5)

%     drawnow
%     yl=ylim*1.1;
%plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
xlim([0,60])
ylim([-1 1])
xlabel('$\theta_{lim}$')
ylabel(['E($\Phi$)'])
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
legend('$\Phi=0$','$\Phi=\pi/2$')
set(gca,'fontsize',font_size_global)

%%
unc_lim = 0.2;
phi_indx_2 = 5;
E_dif = (E_lambda(:,:,1)-E_lambda(:,:,phi_indx_2));
err_mask=sqrt(E_unc(:,:,1).^2+E_unc(:,:,phi_indx_2).^2)<unc_lim;
[valt,indx]=max(E_dif(err_mask))
[xind,yindx]=find(E_dif==valt)

bound_mask = E_dif>sqrt(2);
Violation_significance = (E_dif-sqrt(2))./sqrt(E_unc(:,:,1).^2+E_unc(:,:,phi_indx_2).^2);
[vald,indx]=max(Violation_significance(Violation_significance>0))
[xind,yind]=find(Violation_significance==vald)

%% 
l_indx = 5;
ang_indx = 7;
    x = phi_vec;
%     y = top_corr_bb_vec;%
%     numerator_E = (top_corr_bb_vec+btm_corr_bb_vec-btw_1_corr_bb_vec-btw_2_corr_bb_vec);
%     denominator_E = top_corr_bb_vec+btm_corr_bb_vec+btw_1_corr_bb_vec+btw_2_corr_bb_vec;
%     E_calc = (top_corr_bb_vec+btm_corr_bb_vec-btw_1_corr_bb_vec-btw_2_corr_bb_vec)./(top_corr_bb_vec+btm_corr_bb_vec+btw_1_corr_bb_vec+btw_2_corr_bb_vec);
%     E_raw=(g2_raw.g12.val(:)+g2_raw.g34.val(:)-g2_raw.g14.val(:)-g2_raw.g23.val(:))./(g2_raw.g12.val(:)+g2_raw.g34.val(:)+g2_raw.g14.val(:)+g2_raw.g23.val(:));
    %     y = btm_corr_bb_vec;%
%     if do_g2_err
%         wt = top_corr_bb_unc./2;
%         wb = btm_corr_bb_vec_unc./2;
%         wbt1 = btw_1_corr_bb_vec_unc./2;
%         wbt2 = btw_2_corr_bb_vec_unc./2;
%         err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
%         E_err = abs(E_calc).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
%     else
%         w=y./20;
%         wt = w;
%         wb = w;
%         wbt1 = w;
%         wbt2 = w;
%         err_tot = sqrt(wt.^2+wb.^2+wbt1.^2+wbt2.^2);
%         E_err = abs(E_calc).*sqrt((err_tot./numerator_E).^2+(err_tot./denominator_E).^2);
%     end
    
    
    
    
    xp = linspace(0,max(phi_vec));
    % fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(6)).*(cos(x.*b(5) + 2*pi/b(3))) + b(4);    % Function to fit
%     fit = @(b,x)  b(1).*cos(x + b(2)) + b(3);    % Function to fit
%     best_fit = fitnlm(x,y,fit,[2,0,2],'CoefficientNames',{'Amp','Phase','Offset'}); %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
%     [ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    
    stfig('g2 against phase mask');
    clf
    hold on
%     plot(xp,ysamp_val,'r','LineWidth',1.5)
%     %     drawnow
%     yl=ylim*1.1;
%     plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
    colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
    %     errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    %         'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
    for ii = 1:size(out_corrs,2)
        y12(ii) = out_corrs{ang_indx,ii}.g12.in_shot_corr.(corr_density)(l_indx);
        y34(ii) = out_corrs{ang_indx,ii}.g34.in_shot_corr.(corr_density)(l_indx);
        y23(ii) = out_corrs{ang_indx,ii}.g23.in_shot_corr.(corr_density)(l_indx);
        y14(ii) = out_corrs{ang_indx,ii}.g14.in_shot_corr.(corr_density)(l_indx);
    end
    ht=errorbar(x+0.1,y12,squeeze(wt(l_indx,ang_indx,:)),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hb=errorbar(x,y34,squeeze(wb(l_indx,ang_indx,:)),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    hbt1=errorbar(x+0.1,y23,squeeze(wbt1(l_indx,ang_indx,:)),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);  
    hbt2=errorbar(x,y14,squeeze(wbt2(l_indx,ang_indx,:)),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
    %     legend([hbt1],{'g23'})
    legend([ht hb hbt1 hbt2],{'$C_{12}$', '$C_{34}$','$C_{23}$','$C_{14}$'})
    %     scatter(x,y,'o')
    xlabel('$\Phi$')
    ylabel('$C_{ij}$')
    grid
    box on
    set(gca,'FontSize',19)
    
    stfig('E against phase mask');
%     fit = @(b,x)  b(1).*cos(1.*x + b(2)).^2+1;    % Function to fit
%     best_fit_E = fitnlm(x,E_calc,fit,[1,0.0],'CoefficientNames',{'Amp','Phase'}); %one cos [1.229,1,0.8088,0.906] two cos [1.829,0.01,0.8088,0.906,1.0,0.406]
%     [ysamp_val,ysamp_ci]=predict(best_fit_E,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
    clf
%     hold on
%     plot(xp,ysamp_val,'r','LineWidth',1.5)
    %     drawnow
    %     yl=ylim*1.1;
%     plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
    colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
    %     errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    %         'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
    hE=errorbar(x,squeeze(E_lambda(l_indx,ang_indx,:)),squeeze(E_unc(l_indx,ang_indx,:)),'o','CapSize',0,'MarkerSize',5,'LineWidth',2.5);
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
    
    %     xlim([0,max(phi_vec)])
    %     ylim([1 4])
