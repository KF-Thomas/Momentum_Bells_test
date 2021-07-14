%% data vs phase

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
% opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\full_interferometer\rarity-tapster\tighter_trap\';
%data_folder = '';
log_folder = 'log_Phi.txt';
data_folders = {
    'misc\20210325_k=0,-1_mirror_test_tight_trap\gaussian_test'
    };

force_reimport = true;
force_reimport_norm = true;

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.5e3;%9e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper = Inf;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_norm = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper_norm = 5;%2;%10;%0;% %max allowed number in halo 10

opts.halo_N_lim_both = -1;

y_cut = 11e-3;

% opts.import.shot_num = 1:16; %can select which shots you want to import

%% Calibration settings
L = [-0.4,-0.12].*2;%[-0.3,0.3];%[-0.15,0.15];%[-0.1,0.1];%[-0.618,-0.385];%[-0.5,-0.1];%
L_theta =[-2.*pi, 2.*pi];%[-0.5,0.5];%
norm_folders = [];%[];%
plot_fit = true;
do_g2=true;
do_g2_err=false;
do_range_cut=false;
out_data = {};
phi_vec = [];

num_shots_norm = 0;

%% Run over each folder
for folder_indx = 1:length(data_folders)
    %% import raw data
    data_folder = data_folders{folder_indx};
    opts.import.dir = fullfile(opts.data_root, data_folder);
    opts.logfile = fullfile(opts.import.dir,log_folder);
    if ~ismember(folder_indx,norm_folders)
        %         phi_logs = table2array(readtable(opts.logfile));
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
    slice_type = 'narrow';
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
    
    radius_lim = [0.79,1.17];%[0.61,1.26];%[0.89,1.11];%[0.89,1.16];%[0.9,1.05];%
    
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
    opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
    opts.index.filename = sprintf('index__%s__%.0f', opts.data_dir);
    opts.label = data_folder;
    opts.tag = 0;
    opts.full_out = false;
    opts.bounds = [-0.03, 0.03; -0.03, 0.03];%spacecial bounds
    opts.shot_bounds = [];
    %     if ~exist(opts.fig_dir, 'dir')
    %         mkdir(opts.fig_dir);
    %     end
    % Run the function!
    
    %% Set up out dir
    %set up an output dir %https://gist.github.com/ferryzhou/2269380
    %     if (exist([opts.data_src, '\out'], 'dir') == 0), mkdir(fullfile(opts.data_src, '\out')); end
    %     %make a subfolder with the ISO timestamp for that date
    %     anal_out.dir = sprintf('%sout\\%s\\', ...
    %         [opts.data_src, '\'], datestr(datetime('now'), 'yyyymmddTHHMMSS'));
    %     if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
    
    %% import raw data
    [data, ~] = import_mcp_tdc_data(opts.import);
    
    %% remove any ringing or hotspots
    %     data_ht_spot=hotspot_mask(data);
    %     data.counts_txy=data_ht_spot.masked.counts_txy;
    %     data.num_counts=data_ht_spot.masked.num_counts;
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
    opts.cent.top.threshold = [130,6000,6000].*1e3;
    opts.cent.top.min_threshold = [0,3,3].*1e3;%[16,7,10].*1e3;
    opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.top.method = {'margin','average','average'};
    
    opts.cent.mid.visual = 0; %from 0 to 2
    opts.cent.mid.savefigs = 0;
    opts.cent.mid.threshold = [130,6000,6000].*1e3;
    opts.cent.mid.min_threshold = [0,3,3].*1e3;%[16,7,10].*1e3;
    opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.mid.method = {'margin','average','average'};
    
    opts.cent.btm.visual = 0; %from 0 to 2
    opts.cent.btm.savefigs = 0;
    opts.cent.btm.threshold = [130,6000,6000].*1e3;%[130,2000,2000].*1e3;
    opts.cent.btm.min_threshold = [0,3,3].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
    opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
    opts.cent.btm.method = {'margin','average','average'};
    
    %          opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
%     opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
    opts.cent.t_bounds = {[1.739,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
    %     bec = halo_cent(data_masked,opts.cent);
    num_shots = size(data_masked.shot_num,2);
    
    axis = 1;
    for bec_idx = 1:3
    for idx = 1:num_shots
        lims = [opts.cent.t_bounds{bec_idx}; -0.022, 0.011; -0.08, 0.018];
        this_txy = data_masked.counts_txy{idx};
        trim_txy = masktxy_square(this_txy, lims);
        this_axis = trim_txy(:, axis);
        this_sigma = 1e-5;
        count_hist = smooth_hist(this_axis,'sigma',this_sigma);
        flux = count_hist.count_rate.smooth;
        bin_centres = count_hist.bin.centers;
        gauss =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2));
        kt=1.042e4;%1.607514367076057e+04;%1.038e4;
        gauss_sin =  @(b,x) b(1).*exp(-((x-b(2)).^2)./(2*b(3).^2)).*(1-sin(kt.*x-b(4)))./2;
        fo = statset('TolFun',10^-7,...
            'TolX',1e-5,...
            'MaxIter',1e5,...
            'UseParallel',1);
        %mean(lims(1,:))
        inital_guess_gauss=[max(flux),nanmean(flux.*bin_centres)./mean(flux),1e-3];
        fit_gauss=fitnlm(bin_centres,flux,gauss,inital_guess_gauss,'Options',fo);
        fit_params_gauss = fit_gauss.Coefficients.Estimate';
        
        inital_guess_sin=[0];%[1.037e4,0.0];
        just_sin = @(b,x) gauss_sin([fit_params_gauss,kt,b(1)],x);%gauss_sin([fit_params_gauss,b(1:2)],x);
        fit_sin=fitnlm(bin_centres,flux,just_sin,inital_guess_sin,'Options',fo);
        fit_params_sin = fit_sin.Coefficients.Estimate';
        
        inital_guess_both=[fit_params_gauss,fit_params_sin];
        fit=fitnlm(bin_centres,flux,gauss_sin,inital_guess_both,'Options',fo);
        fit_params = fit.Coefficients.Estimate;
        [prediction,ci]=predict(fit,bin_centres,'Alpha',1-erf(1/sqrt(2)));
        if fit.Rsquared.Adjusted > fit_sin.Rsquared.Adjusted && fit.Rsquared.Adjusted > fit_gauss.Rsquared.Adjusted
            params(idx,bec_idx,:) = fit_params;
        elseif fit_sin.Rsquared.Adjusted > fit.Rsquared.Adjusted
            params(idx,bec_idx,:) = inital_guess_both;
            
        end
    end
    end
    
    %%
    l_t=phi_logs(:,2);
    d_t=data.time_create_write(:,1);
    phi_log_matched = zeros(length(d_t),1).*nan;
    phi_log_check = zeros(length(d_t),1);
    for ii = 1: length(l_t)
        phi_c = phi_logs(ii,3);
        l_c=l_t(ii);
        t_mask=l_c+17<d_t & l_c+29.5>d_t;
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
    is_shot_good = num_check & bec.centre_OK_mid' & phi_log_check';
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
    if ~ismember(folder_indx,norm_folders)
        phi_logs_masked = phi_log_matched(is_shot_good,:);
    end
    %% Find the velocity widths
    opts.bec_width.g0 = const.g0;
    opts.bec_width.fall_time = 0.417;
    
    %% convert data to velocity
    % zero velocity point
    %     t0 = bec_masked_halo.centre_top(:,1);%ones(size(bec_masked_halo.centre_top,1),1).*3.8772;%72;%
    %     x0 = bec_masked_halo.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
    %     y0 = bec_masked_halo.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;
    
    
    
end