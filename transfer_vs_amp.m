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
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '';
data_folders = {
    '20200901_k=0,-1_mirror_vs_amp\Pamp_0'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_0_5'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_1'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_2'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_3'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_4'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_5'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_6'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_7'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_8'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_9'
    
    '20200901_k=0,-1_mirror_vs_amp\Pamp_11'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_12'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_13'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_14'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_15'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_15_t2'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_16'
    '20200901_k=0,-1_mirror_vs_amp\Pamp_19'
    };
% '20200901_k=0,-1_mirror_vs_amp\Pamp_0_25'
% '20200901_k=0,-1_mirror_vs_amp\Pamp_10'

% % Hamming Sinc Pulse with T=6e-6
% % P_amplitude = [1,];
%
% transfer_percentage = [
%     %Pamp, ang, perc
%     1,      0,  0.08333
%     1,      -0.3847,
%     5,
%                     ];

opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10

%% Run over each folder
for this_folder = 1:numel(data_folders)
    %% import raw data
    data_folder = data_folders{this_folder};
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
    opts.cent.correction = 1;
    opts.cent.correction_opts.plots = 0;
    
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
    
    opts.cent.btm.visual = 2; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [30,5000,5000].*1e3;%[50,5000,5000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [0,0,0].*1e3;%[0,3,1.5].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.btm.method = {'margin','average','average'};
    
    opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
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
    
    nbins=50;
    theta_bins = linspace(-pi,pi,nbins);
    phi_bins = linspace(-pi/2,pi/2,nbins);
    
    [theta_top,~] = cart2pol(v_top_zxy(:,2),v_top_zxy(:,3));
    phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
    
    [theta_btm,~] = cart2pol(v_btm_zxy(:,2),v_btm_zxy(:,3));
    phi_btm = atan(v_btm_zxy(:,1)./sqrt(v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2));
    
    for ii = 1:(nbins-1)
        r_btm_zxy_masked = v_btm_zxy(theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1),:);
        r_top_zxy_masked = v_top_zxy(theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1),:);
        v_btm_dens(ii,1) = size(r_btm_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
        v_top_dens(ii,1) = size(r_top_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
        
        r_btm_zxy_masked = v_btm_zxy(phi_bins(ii)<phi_btm & phi_btm<=phi_bins(ii+1),:);
        r_top_zxy_masked = v_top_zxy(phi_bins(ii)<phi_top & phi_top<=phi_bins(ii+1),:);
        v_btm_dens(ii,2) = size(r_btm_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));
        v_top_dens(ii,2) = size(r_top_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));
        
        theta(ii) = mean(theta_bins(ii:(ii+1)));
        phi(ii) = mean(phi_bins(ii:(ii+1)));
    end
    v_btm_dens = v_btm_dens./size(bottom_halo.counts_vel,1);
    v_top_dens = v_top_dens./size(top_halo.counts_vel,1);
    
    %% Append data to structure
    out_data.v_dens.top{this_folder} = v_top_dens;
    out_data.v_dens.btm{this_folder} = v_btm_dens;
    out_data.num_shots(this_folder) = num_shots;
    trans_ratio{this_folder} = v_btm_dens./(v_top_dens+v_btm_dens);
    trans_ratio_unc{this_folder} = v_btm_dens./(v_top_dens+v_btm_dens).*sqrt(1./v_btm_dens+1./(v_top_dens+v_btm_dens));
    equator_trans_ratio(this_folder) = trans_ratio{this_folder}(25,2);
    equator_trans_ratio_unc(this_folder) = trans_ratio_unc{this_folder}(25,2);
    
end
%%
% Eamp = sqrt([0,0.25,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,19]);
Eamp = sqrt([0,0.5,1,2,3,4,5,6,7,8,9,11,12,13,14,15,15,16,19]);

modelfun = @(b,x) b(1).*cos(x(:,1).*b(2)+b(4))+b(3);
fit_1 = fitnlm(Eamp,equator_trans_ratio,modelfun,[-0.5,2*pi/6,0.5,0])

stfig('Transfer efficency against pulse amplitude comp');
clf
Evec = linspace(0,5,1e4);
[ysamp_val,ysamp_ci]=predict(fit_1,Evec','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'
hold on
plot(Evec,ysamp_val,'k','LineWidth',1.5)
drawnow
yl=ylim*1.1;
plot(Evec,ysamp_ci,'color',[1,1,1].*0.5)

curve1 = ysamp_ci(:,1)';
curve2 = ysamp_ci(:,2)';
x1 = Evec;
x2 = [x1, fliplr(x1)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, 'g');
h.FaceColor = [0.31 0.31 0.32].*2;
h.FaceAlpha = 0.5;
parms1 = fit_1.Coefficients.Estimate;
% equator_trans_ratio_unc  = equator_trans_ratio.*0.1;
errorbar(Eamp,equator_trans_ratio,equator_trans_ratio_unc./sqrt(out_data.num_shots),'kx')
% plot(Evec, modelfun(parms1,Evec'))
hold off
xlabel('Transfer Percentage')
ylabel('Pulse Amplitude')