%% Normalised density against phase

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
    %     '20210121_mach_zender_k=0,-1_phi=0'
    %     '20210121_mach_zender_k=0,-1_phi=2pi'
    %     '20210121_mach_zender_k=0,-1_phi=3pi_2'
    %     '20210121_mach_zender_k=0,-1_phi=3pi_2_v2'
    %     '20210121_mach_zender_k=0,-1_phi=5pi_4'
    %     '20210121_mach_zender_k=0,-1_phi=7pi_4'
    %     '20210121_mach_zender_k=0,-1_phi=pi'
    %     '20210121_mach_zender_k=0,-1_phi=pi_2'
    %     '20210121_mach_zender_k=0,-1_phi=pi_4'
    %     '20210121_mach_zender_k=0,-1_phi=pi_v2'
    %     '20210122_mach_zender_k=0,-1_phi=3pi_4'
    '20210204_k=0,-1_norm'
    
    '20210204_phi_0'
    '20210204_phi_2pi'
    '20210204_phi_3pi_2'
    '20210204_phi_3pi_4'
    '20210204_phi_pi'
    '20210204_phi_pi(2)'
    '20210204_phi_pi_2'
    '20210204_phi_pi_4'
    '20210204_phi_5pi_4'
    '20210204_phi_5pi_2'
    '20210204_phi_7pi_4(2)'
    '20210204_phi_3pi'
    '20210204_phi_9pi_4'
    '20210204_phi_mpi_4'
    };

opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10

%% Calibration settings
% phase = [0,2*pi,3*pi/2,3*pi/2,5*pi/4,7*pi/4,pi,pi/2,pi/4,pi,3*pi/4];
norm_folders = [1]; %folders to be used in calibrating out the effects of BEC bleed
phase = [0,2*pi,3*pi/2,3*pi/4,pi,pi,pi/2,pi/4,5*pi/4,5*pi/2,7*pi/4,3*pi,9*pi/4,-pi/4];
L = 0.3;
out_data = [];

%% Run over each folder
for folder_indx = 1:length(data_folders)
    %% import raw data
    data_folder = data_folders{folder_indx};
    opts.import.dir = fullfile(opts.data_root, data_folder);
    
    %% Chose which halo(s) to analyse
    opts.do_top_halo = 1;% analyse the top halo?
    opts.do_btm_halo = 1;% analyse the bottom halo?
    
    %% Chose if you want to look at a narrow or wide slice of the halo
    slice_type = 'wide';
    if strcmp(slice_type,'narrow')
        opts.vel_conv.top.z_mask = [-0.4,0.4];%
        opts.vel_conv.btm.z_mask = [-0.4,0.4];%in units of radius ([-0.68,0.68])
    elseif strcmp(slice_type,'extra wide')
        opts.vel_conv.top.z_mask = [-0.87,0.87];%
        opts.vel_conv.btm.z_mask = [-0.87,0.87];%in units of radius ([-0.68,0.68])
    else
        opts.vel_conv.top.z_mask = [-0.82,0.82];
        opts.vel_conv.btm.z_mask = [-0.82,0.82];%in units of radius ([-0.68,0.68])
    end
    
    %% Import parameters
    tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
    tmp_ylim=[-35e-3, 35e-3];
    tlim=[0,4];
    opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];
    
    opts.num_lim = 2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
    opts.halo_N_lim = 0;%2;%10;%0;% %minimum allowed number in halo 10
    
    opts.plot_dist = false; %do you want to see all the detailed stuff about the halo distributions
    
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
    
    %% set up relevant constants
    hebec_constants
    
    %% find centers
    opts.cent.visual = 0; %from 0 to 2
    opts.cent.savefigs = 0;
    opts.cent.correction = 1;
    opts.cent.correction_opts.plots = 0;
    
    opts.cent.top.visual = 0; %from 0 to 2
    opts.cent.top.savefigs = 0;
    opts.cent.top.threshold = [130,5000,5000].*1e3;
    opts.cent.top.min_threshold = [16,3,3].*1e3;%[16,7,10].*1e3;
    opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.top.method = {'margin','average','average'};
    
    opts.cent.mid.visual = 0; %from 0 to 2
    opts.cent.mid.savefigs = 0;
    opts.cent.mid.threshold = [130,5000,5000].*1e3;
    opts.cent.mid.min_threshold = [16,3,3].*1e3;%[16,7,10].*1e3;
    opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.mid.method = {'margin','average','average'};
    
    opts.cent.btm.visual = 0; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [130,5000,5000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [16,3,3].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
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
    opts.vel_conv.btm.plot_percentage = 0.95;
    opts.vel_conv.btm.title = 'bottom halo';
    opts.vel_conv.btm.const.g0 = const.g0;
    opts.vel_conv.btm.const.fall_distance = const.fall_distance;
    opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
    opts.vel_conv.btm.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
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
    
    v_btm_dens = [];
    v_top_dens = [];
    v_btm_dens_unc = [];
    v_top_dens_unc = [];
    phi_mask_top = (phi_top<0.154& phi_top>-0.154);
    phi_mask_btm = (phi_btm<0.154& phi_btm>-0.154);
    r_btm_zxy_masked=smooth_hist(theta_btm(phi_mask_btm),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
    r_top_zxy_masked=smooth_hist(theta_top(phi_mask_top),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
    v_btm_dens(:,1) = r_btm_zxy_masked.count_rate.smooth;
    v_top_dens(:,1) = r_top_zxy_masked.count_rate.smooth;
    v_btm_dens_unc(:,1) = sqrt(r_btm_zxy_masked.count_rate.smooth).*sqrt(abs(r_btm_zxy_masked.bin.edge(1:end-1)...
        -r_btm_zxy_masked.bin.edge(2:end)));
    v_top_dens_unc(:,1) = sqrt(r_top_zxy_masked.count_rate.smooth).*sqrt(abs(r_top_zxy_masked.bin.edge(1:end-1)...
        -r_top_zxy_masked.bin.edge(2:end)));
    
    
    r_btm_zxy_masked=smooth_hist(phi_btm,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
    r_top_zxy_masked=smooth_hist(phi_top,'sigma',0.04,'lims',[-pi/2,pi/2],'bin_num',nbins);
    v_btm_dens(:,2) = r_btm_zxy_masked.count_rate.smooth;
    v_top_dens(:,2) = r_top_zxy_masked.count_rate.smooth;
    
    v_top_dens_2d = hist3([theta_top phi_top],'Nbins',[nbins nbins]);
    v_btm_dens_2d = hist3([theta_btm phi_btm],'Nbins',[nbins nbins]);
    
    theta = linspace(-pi,pi,nbins);
    phi = linspace(-pi/2,pi/2,nbins);
    
    %% ratio in spherical coordinates
    top_dens_norm = (v_top_dens(:,2))./(v_btm_dens(:,2)+v_top_dens(:,2));
    if ismember(folder_indx,norm_folders)
        if ~exist('cal_dens_top','var')
            cal_dens_top = zeros(nbins,1);
            cal_dens_btm = zeros(nbins,1);
        end
        cal_dens_top = cal_dens_top + v_top_dens(:,2)./size(norm_folders,1);
        cal_dens_btm = cal_dens_btm + v_btm_dens(:,2)./size(norm_folders,1);
    else
        phi_mask = phi<L/2 & phi>-L/2;
        top_dens_avg = trapz(phi(phi_mask),top_dens_norm(phi_mask))./range(phi(phi_mask));%./L;
        top_dens_std = sqrt(trapz(phi(phi_mask),(top_dens_avg-top_dens_norm(phi_mask)).^2)./range(phi(phi_mask)));
        
        trans_ratio_top = (v_top_dens(:,2)-cal_dens_top)./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_top.*2);
        trans_ratio_btm = (v_btm_dens(:,2)-cal_dens_btm)./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_btm.*2);
        
        trans_ratio_top_avg = trapz(phi(phi_mask),trans_ratio_top(phi_mask))./range(phi(phi_mask));
        trans_ratio_top_unc = sqrt(trapz(phi(phi_mask),(trans_ratio_top_avg-trans_ratio_top(phi_mask)).^2)./range(phi(phi_mask)));
        
        trans_ratio_btm_avg = trapz(phi(phi_mask),trans_ratio_btm(phi_mask))./range(phi(phi_mask));
        trans_ratio_btm_unc = sqrt(trapz(phi(phi_mask),(trans_ratio_btm_avg-trans_ratio_btm(phi_mask)).^2)./range(phi(phi_mask)));
        
        out_data(folder_indx-size(norm_folders,1),:) = [top_dens_avg,top_dens_std];
        out_data_trans_btm(folder_indx-size(norm_folders,1),:) = [trans_ratio_btm_avg,trans_ratio_btm_unc];
        out_data_trans_top(folder_indx-size(norm_folders,1),:) = [trans_ratio_top_avg,trans_ratio_top_unc];
        out_data_vec(folder_indx-size(norm_folders,1),:) = top_dens_norm;
    end
end
%%
% [x, y, w] = combine_data(phase,out_data(:,1),out_data(:,2));
[x, y, w] = combine_data(phase,1-out_data_trans_btm(:,1),out_data_trans_btm(:,2));
% [x, y, w] = combine_data(phase,out_data_trans_top(:,1),out_data_trans_top(:,2));
[x_vec, y_vec] = combine_data(phase,out_data_vec);

stfig('dens over sphere');
[a, b] = sort(x_vec);
pcolor(a,phi,y_vec(b,:)')
% pcolor(phi,x_vec,y_vec)
% caxis([0.28 0.6])
caxis([0.3 0.7])
ylim([-0.3,0.3])
ylabel('$\phi$')
xlabel('phase')
colorbar

xp = linspace(min(x),1.1*max(x));
fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(3)) + b(4);    % Function to fit
best_fit = fitnlm(x,y,fit,[0.329,1,-0.8088,0.5],'CoefficientNames',{'Amp','Freq','Phase','Offset'});
[ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'

stfig('normalised density of top halo against phase');
clf
hold on
plot(xp,ysamp_val,'r','LineWidth',1.5)
drawnow
yl=ylim*1.1;
plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]./255;
errorbar(x,y,w,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
xlabel('$\phi$')
ylabel('normalised density top halo')
grid
box on
set(gca,'FontSize',19)
xlim([xp(1),xp(end)])