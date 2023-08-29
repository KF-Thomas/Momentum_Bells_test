%% Transfer over amplitude

%% Initializing path
clear all;
% close all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

%% Import directories
% opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\pulse_characterisation\';
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '20230406,k=0,+1, p=0.3,T=var,f=18_delay=210_AB';%'20230406,k=0,+1, p=0.3,T=40,f=18_delay=210_AB';
log_folder = 'log_Phi.txt';
log_lab_folder = 'log_LabviewMatlab.txt';
data_folders = {
% ''
    %original scan
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_0_t2'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_0_25_t2'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_0_5'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_1'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_2'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_3'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_4'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_5'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_6'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_7'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_8'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_9'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_10_t2'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_11'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_12'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_13'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_14'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_15'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_15_t2'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_16'
%     '20200901_k=0,-1_transfer_vs_amp\Pamp_19_t2'
%     
%     '20200907_k=-1,-2_transfer_vs_amp\Pamp_0'
%     '20200907_k=-1,-2_transfer_vs_amp\Pamp_5_3'
%     '20200907_k=-1,-2_transfer_vs_amp\Pamp_11_5'

    %new scan
%     '20210201_bragg_pulse_testing\scan\amp_0'
%     '20210201_bragg_pulse_testing\scan\amp_2'
%     '20210201_bragg_pulse_testing\scan\amp_2_5'
%     '20210201_bragg_pulse_testing\scan\amp_3\alpha_3_5'
%     '20210201_bragg_pulse_testing\scan\amp_3_5\alpha_3_5'
%     '20210201_bragg_pulse_testing\scan\amp_4'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_3_5'
%     '20210201_bragg_pulse_testing\scan\amp_5'
%     '20210201_bragg_pulse_testing\scan\amp_5_5'
    
    %alpha scan
%     '20210201_bragg_pulse_testing\scan\amp_3_5\alpha_3'
%     '20210201_bragg_pulse_testing\scan\amp_3_5\alpha_3_5'
%     '20210201_bragg_pulse_testing\scan\amp_3_5\alpha_4'
%     '20210201_bragg_pulse_testing\scan\amp_3_5\alpha_4_5'

%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_2_5'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_3'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_3_25'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_3_5'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_3_75'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_4'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_4_5'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_5'
%     '20210201_bragg_pulse_testing\scan\amp_4_5\alpha_6'
    };

%
% '20200901_k=0,-1_mirror_vs_amp\Pamp_10'
% '20200901_k=0,-1_mirror_vs_amp\Pamp_0_25'

% % Hamming Sinc Pulse with T=6e-6
% % P_amplitude = [1,];
%
% transfer_percentage = [
%     %Pamp, ang, perc
%     1,      0,  0.08333
%     1,      -0.3847,
%     5,
%                     ];

opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.0e3;%2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10

%% Calibration settings
% norm_folders = [1 22]; %folders to be used in calibrating out the effects of BEC bleed
norm_folders = [1]; %folders to be used in calibrating out the effects of BEC bleed
% btm_halo_comp = [22 23 24];%transfering the bottom halo for comparison
btm_halo_comp = [];
% Eamp = sqrt([0,0.25,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,19,0,5.3,11.5]); %amplitude vector %sqrt([0,5.3,11.5]);%
Eamp = 1:6;%[0,0.23,0.26,0.28,0.29];%[0,0.085,0.115];%[0,0.1,0.145,0.18];%[0,0.23,0.26,0.28,0.29];%[0,linspace(0.1,0.26,9-1)];%[0,0.232];%[0,0.085,0.14,0.15,0.155];%[0,0.07,0.14,0.15];%[0,linspace(0.1,0.29,12-1)];%[0,1];%[0,linspace(0.16,0.28,12-1)];%linspace(0,0.6,12);%[0,linspace(0.16,0.28,4)];%[0,0.52,0.55];%[0,linspace(0.35,0.65,5-1)];%[0,linspace(0.35,0.65,10-1)];%linspace(0,20,10).*1e3;%[0,2,2.5,3,3.5,4,4.5,5,5.5];
% Eamp = sqrt([0,3,3.5,4,4.5]);
% Eamp = sqrt([0,2.5,3,3.25,3.5,3.75,4,4.5,5,6]);
num_norms = length(norm_folders);
folder_indxs = [norm_folders, 1:numel(data_folders)];
calibration_density_btm= [];
calibration_shot_num = [];
cal_dens_btm = 0;

opts.logfile = fullfile(fullfile(opts.data_root, data_folder),log_folder);
    opts.loglabfile = fullfile(fullfile(opts.data_root, data_folder),log_lab_folder);
    log_check = isfile(opts.logfile);
if log_check
    phi_logs = table2array(readtable(opts.logfile));
    unique_phi = unique(phi_logs(:,3));
    folder_indxs = [norm_folders, 1:numel(unique_phi)];
end

loop_num = length(folder_indxs);

%% Run over each folder
cal_dens_top_data_top = zeros(151,2);
cal_dens_top_data_btm = zeros(151,2);
for folder_indx = 1:loop_num
    
    %% import raw data

    this_folder = folder_indxs(folder_indx);

    if log_check
        phi_mask = phi_logs(:,3) == unique_phi(this_folder);
        opts.import.shot_num = phi_logs(phi_mask,1).'; %
    else
        data_folder = data_folders{this_folder};
        opts.import = rmfield(opts.import,'shot_num');
    end

    
    opts.import.dir = fullfile(opts.data_root, data_folder);
    [data, ~] = import_mcp_tdc_data(opts.import);
    
    %% remove any ringing
    opts.ring_lim = 0.09e-6;%0.1e-6;%-1;%0;%0.101 %how close can points be in time
    data_masked = ring_removal(data,opts.ring_lim);
    
    %% set up relevant constants
    hebec_constants
    
    %% find centers
    opts.cent.visual = 0; %from 0 to 2
    opts.cent.savefigs = 0;
    opts.cent.correction = 2;
    opts.cent.correction_opts.plots = 0;
    opts.cent.nan_cull = 1;
    
    opts.cent.top.visual = 0; %from 0 to 2
    opts.cent.top.savefigs = 0;
    opts.cent.top.threshold = [130,2000,2000].*1e3;
    opts.cent.top.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.top.method = {'margin','average','average'};
    
    opts.cent.mid.visual = 0; %from 0 to 2
    opts.cent.mid.savefigs = 0;
    opts.cent.mid.threshold = [140,2000,2000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.mid.min_threshold = [16,9,9].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.mid.method = {'margin','average','average'};
    
    opts.cent.btm.visual = 0; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [30,5000,5000].*1e3;%[50,5000,5000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [0,0,0].*1e3;%[0,3,1.5].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.btm.method = {'margin','average','average'};
    
    opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.884,3.896],[3.75,4]};%
%     opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
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
    is_shot_good = is_shot_good & bec.centre_OK_top';
    is_shot_good = is_shot_good & bec.centre_OK_btm';
    num_shots = sum(is_shot_good);
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
    opts.vel_conv.top.v_mask=[0.89,1.11].*0.065; %bounds on radisu as multiple of radius value
    opts.vel_conv.top.z_mask = [-0.8,0.8]; %in units of radius (standard [-0.76,0.76])
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
    opts.vel_conv.btm.v_mask=[0.89,1.11].*0.065; %bounds on radisu as multiple of radius value
    opts.vel_conv.btm.z_mask = [-0.8,0.8]; %in units of radius
    opts.vel_conv.btm.y_mask = [-1.9,1.9]; %in units of radius
    opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point
    
    opts.vel_conv.btm.centering_correction = [0,0,0]; %correctoin shift to the centering in m/s
    
    opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
    opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
    opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
    opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;
    
    %%
    bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);
    
    %% Mask out halos with nums to low
    halo_N_check_top = top_halo_intial.num_counts>opts.halo_N_lim;
    halo_N_check_btm = bottom_halo_intial.num_counts>opts.halo_N_lim;
    halo_N_check = halo_N_check_top & halo_N_check_btm;
    top_halo = struct_mask(top_halo_intial,halo_N_check);
    bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
    bec_halo = struct_mask(bec_masked_halo,halo_N_check);
    
    %%
    v_top_zxy = cell2mat(top_halo.counts_vel);
    v_btm_zxy = cell2mat(bottom_halo.counts_vel);
    
    %     nbins=50;%50
    %     theta_bins = linspace(-pi,pi,nbins);
    %     phi_bins = linspace(-pi/2,pi/2,nbins);
    
    [theta_top,~] = cart2pol(v_top_zxy(:,2),v_top_zxy(:,3));
    phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
    
    [theta_btm,~] = cart2pol(v_btm_zxy(:,2),v_btm_zxy(:,3));
    phi_btm = atan(v_btm_zxy(:,1)./sqrt(v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2));
    
    r_btm_zxy_masked=smooth_hist(theta_btm,'sigma',0.1,'lims',[-pi,pi],'bin_num',151);
    r_top_zxy_masked=smooth_hist(theta_top,'sigma',0.1,'lims',[-pi,pi],'bin_num',151);
    v_btm_dens(:,1) = r_btm_zxy_masked.count_rate.smooth;
    v_top_dens(:,1) = r_top_zxy_masked.count_rate.smooth;
    
    
    r_btm_zxy_masked=smooth_hist(phi_btm,'sigma',0.1,'lims',[-pi/2,pi/2],'bin_num',151);
    r_top_zxy_masked=smooth_hist(phi_top,'sigma',0.1,'lims',[-pi/2,pi/2],'bin_num',151);
    v_btm_dens(:,2) = r_btm_zxy_masked.count_rate.smooth;
    v_top_dens(:,2) = r_top_zxy_masked.count_rate.smooth;
    


    if folder_indx == num_norms+1 && num_norms>0
        cal_dens_top_data_top = squeeze(nansum(calibration_density_top.*calibration_shot_num,1)./nansum(calibration_shot_num(:,1,1)));
        cal_dens_top_data_btm = squeeze(nansum(calibration_density_btm.*calibration_shot_num,1)./nansum(calibration_shot_num(:,1,1)));
        
        
%         cal_dens_btm_data_top = squeeze(nansum(calibration_density_btm_data_top.*calibration_shot_num_btm_data,1)./nansum(calibration_shot_num_btm_data(:,1,1)));
%         cal_dens_btm_data_btm = squeeze(nansum(calibration_density_btm_data_btm.*calibration_shot_num_btm_data,1)./nansum(calibration_shot_num_btm_data(:,1,1)));
        
    end
    v_btm_dens = v_btm_dens./size(bottom_halo.counts_vel,1);
    %     stfig('dens btm comp')
    %     hold on
    %     plot(phi,v_btm_dens(:,2))
    %     modelfun = @(b,x) b(1).*exp(-b(2).*(sin(x)-b(3)).^2);
    %     fit_0 = fitnlm(phi(1:35),v_btm_dens(1:35,2),modelfun,[100,0.5,1])
    %     parms0 = fit_0.Coefficients.Estimate;
    %     plot(phi, modelfun([parms0],phi),'--','LineWidth',2.0)
    %     hold off
    v_top_dens = v_top_dens./size(top_halo.counts_vel,1);
    
    
    %% Append data to structure
    if folder_indx>num_norms
        out_data.v_dens.top{this_folder} = v_top_dens;
        out_data.v_dens.btm{this_folder} = v_btm_dens;
        out_data.num_shots(this_folder) = num_shots;
        out_data.num_counts.top(this_folder) = size(v_top_zxy,1);
        out_data.num_counts.btm(this_folder) = size(v_btm_zxy,1);
        
        N_Bec(this_folder) = nanmean(bec_halo.counts_mid);
        Cen_test(this_folder,:) = nanmean(bec_halo.width_mid);
        
        if ismember(this_folder,btm_halo_comp)
            cal_dens_top = cal_dens_btm_data_top;
            cal_dens_btm = cal_dens_btm_data_btm;
            N_cal = N_Bec(22);
            cal_N = cal_N_btm_data;
        else
            cal_dens_top = cal_dens_top_data_top;
            cal_dens_btm = cal_dens_top_data_btm;
            N_cal = N_Bec(1);
%             cal_N = cal_N_top_data;
        end
        
%         if Eamp(this_folder).^2<3
%             Cbtm = N_Bec(this_folder)./N_cal;
%             Ctop = (cal_N./((size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots)).^(-1);
%         elseif Eamp(this_folder).^2<9
%             Cbtm = N_Bec(this_folder)./N_cal;
%             Ctop = (cal_N./((size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots)).^(-1);
%         else
%             Cbtm = (cal_N./((size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots)).^(-1);
%             Ctop = (cal_N./((size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots)).^(-1);
%         end
        %         Ctop = Cbtm;
        
        
        %(v_top_dens+v_btm_dens)./(cal_dens_top+cal_dens_btm); %sum
        %normalisation
        
        %(cal_N./((size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots)).^(-1);
        %%number normalisation
        
        %Bec number normalisation
        %N_test(this_folder)./N_test(1)
        
        %Bec width normalisation
        %(Cen_test(1,1)./Cen_test(this_folder,1)).^1.5;
        
        %1;%no normalisation (unadjusted calibration)
        Ctop = 1;
        Cbtm = 1;
        C_vec(this_folder,:) = [Ctop,Cbtm];
        trans_ratio_top{this_folder} = (v_top_dens-Ctop.*cal_dens_top)./(v_top_dens+v_btm_dens-Ctop.*cal_dens_top.*2);
        trans_ratio_btm{this_folder} = (v_btm_dens-Cbtm.*cal_dens_btm)./(v_top_dens+v_btm_dens-Cbtm.*cal_dens_btm.*2);
        trans_ratio_naive{this_folder} = (v_btm_dens)./(v_top_dens+v_btm_dens);%
        
        nbins = size(v_btm_dens,1);
        trans_ratio_unc{this_folder} = v_btm_dens./(v_top_dens+v_btm_dens).*sqrt(1./v_btm_dens+1./(v_top_dens+v_btm_dens));
%         equator_trans_ratio_top(this_folder) = trans_ratio_top{this_folder}(round(nbins/2),2);
%         equator_trans_ratio_btm(this_folder) = trans_ratio_btm{this_folder}(round(nbins/2),2);
%         equator_trans_naive(this_folder) = trans_ratio_naive{this_folder}(round(nbins/2),2);
%         equator_trans_ratio_unc(this_folder) = trans_ratio_unc{this_folder}(round(nbins/2),2);
        indx_vec = 68:84;%65:87;%:73;
        equator_trans_ratio_top(this_folder) = nanmean(trans_ratio_top{this_folder}(indx_vec,2));
        equator_trans_ratio_btm(this_folder) = nanmean(trans_ratio_btm{this_folder}(indx_vec,2));
        equator_trans_naive(this_folder) = nanmean(trans_ratio_naive{this_folder}(indx_vec,2));
        equator_trans_ratio_unc(this_folder) = nanmean(trans_ratio_unc{this_folder}(indx_vec,2));

        equator_trans_ratio_btm_unc(this_folder) = nanstd(trans_ratio_btm{this_folder}(indx_vec,2));
        equator_trans_ratio_top_unc(this_folder) = nanstd(trans_ratio_top{this_folder}(indx_vec,2));
        equator_trans_naive_unc(this_folder) = nanstd(trans_ratio_naive{this_folder}(indx_vec,2));
    else
        if ismember(this_folder,btm_halo_comp)
            calibration_density_btm_data_btm(folder_indx,:,:) = v_btm_dens;
            calibration_density_btm_data_top(folder_indx,:,:) = v_top_dens;
            calibration_shot_num_btm_data(folder_indx,:,:) = size(bottom_halo.counts_vel,1).*ones(size(v_btm_dens));
            cal_N_btm_data = (size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots;
            
        else
            calibration_density_btm(folder_indx,:,:) = v_btm_dens;
            calibration_density_top(folder_indx,:,:) = v_top_dens;
            calibration_shot_num(folder_indx,:,:) = size(bottom_halo.counts_vel,1).*ones(size(v_btm_dens));
            cal_N_top_data = (size(v_top_zxy,1)+size(v_btm_zxy,1))./num_shots;
        end
    end
    
end
%%
Eamp = unique_phi;
% top_indx = setdiff(1:numel(data_folders),btm_halo_comp);
top_indx = setdiff(1:numel(unique_phi),btm_halo_comp);
% top_indx = 1:10;%
% Eamp = sqrt([0,0.25,0.5,1,2,3,4,5,6,7,8,9,11,12,13,14,15,15,16,19]);
% Eamp = sqrt([0,0.25,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,19,5.3]);
equator_trans_ratio_top_bounded = equator_trans_ratio_top;
% equator_trans_ratio_top_bounded>1 |
% equator_trans_ratio_top_bounded(equator_trans_ratio_top_bounded>1 | equator_trans_ratio_top_bounded<0) = nan;
equator_trans_ratio_btm_bounded = equator_trans_ratio_btm;
% equator_trans_ratio_btm_bounded(equator_trans_ratio_btm_bounded>1 | equator_trans_ratio_btm_bounded<0) = nan;

equator_trans_ratio_both = nanmean([equator_trans_ratio_top_bounded(top_indx);equator_trans_ratio_btm_bounded(top_indx)]);
equator_trans_ratio_both_btm = nanmean([equator_trans_ratio_top_bounded(btm_halo_comp);equator_trans_ratio_btm_bounded(btm_halo_comp)]);

%thresholding
equator_trans_ratio_both(equator_trans_ratio_both<0) = 0;
equator_trans_ratio_both(equator_trans_ratio_both>1) = 1;

% trans_naive = [out_data.v_dens.btm{:}]./([out_data.v_dens.top{:}]+[out_data.v_dens.btm{:}]);

colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;

% modelfun = @(b,x) b(1).*cos(x(:,1).*b(2)+b(4))+b(3);
% fit_1 = fitnlm(Eamp,equator_trans_ratio,modelfun,[-0.5,2*pi/6,0.5,0])
% Eamp = Eamp.';
modelfun = @(b,x) min(10e6.*ones(size(x)),real(b(1).*sin(b(4)+(b(3).*x(:,1)).^b(2)).^2));
fit_top = fitnlm(Eamp,equator_trans_ratio_top,modelfun,[1,0.5,30,0])
fit_btm = fitnlm(Eamp,equator_trans_ratio_btm,modelfun,[1,0.5,30,pi/2])
% fit_naive = fitnlm(Eamp,equator_trans_naive,modelfun,[1,0,1,1])
% fit_both = fitnlm(Eamp(top_indx),equator_trans_ratio_both,modelfun,[1,0,1,1])

stfig('Transfer efficency against pulse amplitude comp');
clf
grid on
hold on
Evec = linspace(0,0.5,1e4);
[ysamp_val,ysamp_ci]=predict(fit_top,Evec','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
plot(Evec,ysamp_val,'k','LineWidth',1.5)
drawnow
yl=ylim*1.1;
% plot(Evec,ysamp_ci,'color',[1,1,1].*0.5)

% curve1 = ysamp_ci(:,1)';
% curve2 = ysamp_ci(:,2)';
% x1 = Evec.^2;
% x2 = [x1, fliplr(x1)];
% inBetween = [curve1, fliplr(curve2)];
% h = fill(x2, inBetween, 'g');
% h.FaceColor = [0.31 0.31 0.32].*2;
% h.FaceAlpha = 0.5;
% parms1 = fit_top.Coefficients.Estimate;
% equator_trans_ratio_unc  = equator_trans_ratio.*0.1;
% p_top = errorbar(Eamp(top_indx).^2,equator_trans_ratio_both,0.5./sqrt(out_data.num_shots(top_indx)),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
%     'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
% % 
% p_btm = errorbar(Eamp(btm_halo_comp).^2,equator_trans_ratio_both_btm,0.5./sqrt(out_data.num_shots(btm_halo_comp)),'o','CapSize',0,'MarkerSize',5,'Color',[0.7,0.4,.3],...
%     'MarkerFaceColor',[0.8,0.2,.2],'LineWidth',2.5)

p_nav = errorbar(Eamp,equator_trans_naive,equator_trans_naive_unc,'o','CapSize',0,'MarkerSize',5,'Color',[0.5,0.5,.2],...
    'MarkerFaceColor',[0.5,0.7,.05],'LineWidth',2.5);
% %
p_btm = errorbar(Eamp,equator_trans_ratio_btm,equator_trans_ratio_btm_unc,'o','CapSize',0,'MarkerSize',5,'Color',[0.8,0.3,.1],...
    'MarkerFaceColor',[0.8,0.2,.2],'LineWidth',2.5);

p_top = errorbar(Eamp,equator_trans_ratio_top,equator_trans_ratio_top_unc,'o','CapSize',0,'MarkerSize',5,'Color',[0.5,0.8,.7],...
    'MarkerFaceColor',[0.6,0.8,.7],'LineWidth',2.5);
% % 
% errorbar(Eamp(top_indx).^2,equator_trans_ratio_both,0.5./sqrt(out_data.num_shots(top_indx)),'o','CapSize',0,'MarkerSize',5,'Color',[0.5,0.8,.7],...
%     'MarkerFaceColor',[0.6,0.8,.7],'LineWidth',2.5)


% plot(Evec.^2, modelfun([1,pi/10,2],Evec'))
legend([p_top,p_btm p_nav],'Top Halo transfer','Bottom Halo transfer','Naive Transfer')
set(gca,'FontSize',22)
hold off
ylabel('Transfer Percentage')
xlabel('Pulse Amplitude')
ylim([0 1])

%%
stfig('transfer plot');
clf
plt_indx = 3;
unique_phi(plt_indx)
subplot(1,2,1)
plot(r_top_zxy_masked.bin.centers,out_data.v_dens.btm{plt_indx}(:,2)-cal_dens_btm(:,2),'LineWidth',2.5)
hold on
plot(r_top_zxy_masked.bin.centers,out_data.v_dens.top{plt_indx}(:,2)-cal_dens_top(:,2),'LineWidth',2.5)
% plot(r_top_zxy_masked.bin.centers,out_data.v_dens.btm{plt_indx}(:,2),'LineWidth',2.5)
% hold on
% plot(r_top_zxy_masked.bin.centers,out_data.v_dens.top{plt_indx}(:,2),'LineWidth',2.5)

plot(r_top_zxy_masked.bin.centers(indx_vec),out_data.v_dens.btm{plt_indx}(indx_vec,2)-cal_dens_btm(indx_vec,2),'k--','LineWidth',2.0)
plot(r_top_zxy_masked.bin.centers(indx_vec),out_data.v_dens.top{plt_indx}(indx_vec,2)-cal_dens_top(indx_vec,2),'k--','LineWidth',2.0)

ylabel('density')
xlabel('elevation  (rad)')
set(gca,'FontSize',18)
xlim([-0.8, 0.8])
legend('Density of Bottom','Density of Top','Averaged Region')
subplot(1,2,2)
plot(r_top_zxy_masked.bin.centers,trans_ratio_btm{plt_indx}(:,2),'LineWidth',2.5)
hold on
plot(r_top_zxy_masked.bin.centers,trans_ratio_naive{plt_indx}(:,2),'LineWidth',2.5)

plot(r_top_zxy_masked.bin.centers(indx_vec),trans_ratio_btm{plt_indx}(indx_vec,2),'k--','LineWidth',2.0)
plot(r_top_zxy_masked.bin.centers(indx_vec),trans_ratio_naive{plt_indx}(indx_vec,2),'k--','LineWidth',2.0)

xlim([-0.8, 0.8])
ylabel('transfer percentage')
xlabel('elevation (rad)')
set(gca,'FontSize',18)
legend('Transfer from Bottom','Transfer from Top','Averaged Region')
ylim([0 1])

%%
plt_indx = 4;
stfig('flip comp');
clf
plot(r_top_zxy_masked.bin.centers,trans_ratio_btm{plt_indx}(:,2),'LineWidth',2.5)
hold on
plot(r_top_zxy_masked.bin.centers,trans_ratio_naive{plt_indx}(:,2),'LineWidth',2.5)

plot(r_top_zxy_masked.bin.centers,flip(trans_ratio_btm{plt_indx}(:,2)),'--','LineWidth',2.5)
plot(r_top_zxy_masked.bin.centers,flip(trans_ratio_naive{plt_indx}(:,2)),'--','LineWidth',2.5)
xlim([-0.8, 0.8])
xlabel('elevation (rad)')
ylabel('transfer percentage')
set(gca,'FontSize',18)
legend('Transfer from Bottom','Transfer from Top','Transfer from Bottom flipped','Transfer from Top flipped')
ylim([0 1])
%% Plot num checks

stfig('num checks');
clf
subplot(2,1,2)
plot(Eamp,(out_data.num_counts.top(:)+out_data.num_counts.btm(:)),'o')
xlabel('Amp')
ylabel('total counts in both halos')
subplot(2,1,1)
plot(Eamp,(out_data.num_counts.top(:)+out_data.num_counts.btm(:))./out_data.num_shots','o')
xlabel('Amp')
ylabel('avg counts in both halos')
set(gca,'FontSize',17)