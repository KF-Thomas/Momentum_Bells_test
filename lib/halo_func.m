function out = halo_func(data_folder,opts)
%% Background stuff
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
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
opts.bounds = [-0.038, 0.038; -0.038, 0.038];%spacecial bounds
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
    opts.import.dir = dirs_list;
        opts.import.cache_save_dir = fullfile(dirs_list, 'cache', 'import\');

    [data, ~] = import_mcp_tdc_data(opts.import);
end
% [data, ~] = import_mcp_tdc_data(opts.import);

%% remove any ringing
data_ht_spot=hotspot_mask(data);
data.counts_txy=data_ht_spot.masked.counts_txy;
data.num_counts=data_ht_spot.masked.num_counts;
opts.ring_lim = -1;%0.09e-6;%0.1e-6;%0;%0.101 %how close can points be in time
data_masked = ring_removal(data,opts.ring_lim);

%% set up relevant constants
hebec_constants
%%
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
opts.vel_conv.top.v_mask=opts.radius_lim; %bounds on radisu in units of m/s
opts.vel_conv.top.ang_lim = opts.ang_lim; %angular limits of the azimuthal angle
opts.vel_conv.top.vxy_mask = [1i,0.85];
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
opts.vel_conv.btm.v_mask=opts.radius_lim;%[0.61,1.26];%[0.3,1.61];%[0.61,1.26];%[0.89,1.11]; %bounds on radisu as multiple of radius value
opts.vel_conv.btm.ang_lim = opts.ang_lim; %angular limits of the azimuthal angle
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

out.v_top_zxy = cell2mat(top_halo.counts_vel_norm);
out.v_top_zxy_unnorm = cell2mat(top_halo.counts_vel);
out.r_dist_top = sqrt(out.v_top_zxy(:,1).^2+out.v_top_zxy(:,2).^2+out.v_top_zxy(:,3).^2);
out.r_dist_top_unnorm = sqrt(out.v_top_zxy_unnorm(:,1).^2+out.v_top_zxy_unnorm(:,2).^2+out.v_top_zxy_unnorm(:,3).^2);
N_top = top_halo.num_counts;

out.v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm);
out.v_btm_zxy_unnorm = cell2mat(bottom_halo.counts_vel);
if ~isempty(out.v_btm_zxy)
    out.r_dist_btm = sqrt(out.v_btm_zxy(:,1).^2+out.v_btm_zxy(:,2).^2+out.v_btm_zxy(:,3).^2);
    out.r_dist_btm_unnorm = sqrt(out.v_btm_zxy_unnorm(:,1).^2+out.v_btm_zxy_unnorm(:,2).^2+out.v_btm_zxy_unnorm(:,3).^2);
else
    out.r_dist_btm = [];
    out.r_dist_btm_unnorm = [];
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
    [theta_top,out.rxy_top] = cart2pol(out.v_top_zxy_unnorm(:,2),out.v_top_zxy_unnorm(:,3));
    phi_top = atan(out.v_top_zxy(:,1)./sqrt(out.v_top_zxy(:,2).^2+out.v_top_zxy(:,3).^2));
else
    theta_top = [];
    phi_top = [];
end

if opts.do_btm_halo
    [theta_btm,out.rxy_btm] = cart2pol(out.v_btm_zxy_unnorm(:,2),out.v_btm_zxy_unnorm(:,3));
    phi_btm = atan(out.v_btm_zxy(:,1)./sqrt(out.v_btm_zxy(:,2).^2+out.v_btm_zxy(:,3).^2));
else
    theta_btm = [];
    phi_btm = [];
end
for ii = 1:(nbins-1)
    r_btm_zxy_masked = out.r_dist_btm_unnorm(theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1));
    r_top_zxy_masked = out.r_dist_top_unnorm(theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1));
    v_btm_r(ii,1) = mean(r_btm_zxy_masked);
    v_top_r(ii,1) = mean(r_top_zxy_masked);
    %     v_btm_dens(ii,1) = size(r_btm_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
    %     v_top_dens(ii,1) = size(r_top_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
    
    r_btm_zxy_masked = out.r_dist_btm_unnorm(phi_bins(ii)<phi_btm & phi_btm<=phi_bins(ii+1));
    r_top_zxy_masked = out.r_dist_top_unnorm(phi_bins(ii)<phi_top & phi_top<=phi_bins(ii+1));
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

if opts.do_top_halo
    out.r_top_zxy_masked=smooth_hist(theta_top(phi_mask_top),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
    out.v_top_dens(:,1) = out.r_top_zxy_masked.count_rate.smooth./num_shots;
    out.v_top_dens_unc(:,1) = sqrt(out.r_top_zxy_masked.count_rate.smooth).*sqrt(abs(out.r_top_zxy_masked.bin.edge(1:end-1)...
        -out.r_top_zxy_masked.bin.edge(2:end)));
    out.r_top_zxy_masked=smooth_hist(phi_top,'sigma',phi_sig,'lims',[-pi/2,pi/2],'bin_num',nbins);
    out.v_top_dens(:,2) = out.r_top_zxy_masked.count_rate.smooth./num_shots;
    out.top_phi = out.r_top_zxy_masked.bin.centers;
    out.v_top_dens_2d = hist3([theta_top phi_top],'Nbins',[nbins nbins]);

    out.v_hist = smooth_hist(out.v_top_zxy_unnorm(:,1),'sigma',0.2e-2,'lims',[-6.5e-2.*7,6.5e-2.*3],'bin_num',nbins*3);
    out.v_dens = out.v_hist.count_rate.smooth./num_shots;
    out.v_dens_cen =  out.v_hist.bin.centers;
end
if opts.do_btm_halo
    out.r_btm_zxy_masked=smooth_hist(theta_btm(phi_mask_btm),'sigma',0.04,'lims',[-pi,pi],'bin_num',nbins);
    out.v_btm_dens(:,1) = out.r_btm_zxy_masked.count_rate.smooth./num_shots;
    out.v_btm_dens_unc(:,1) = sqrt(out.r_btm_zxy_masked.count_rate.smooth).*sqrt(abs(out.r_btm_zxy_masked.bin.edge(1:end-1)...
        -out.r_btm_zxy_masked.bin.edge(2:end)));
    out.r_btm_zxy_masked=smooth_hist(phi_btm,'sigma',phi_sig,'lims',[-pi/2,pi/2],'bin_num',nbins);
    out.v_btm_dens(:,2) = out.r_btm_zxy_masked.count_rate.smooth./num_shots;
    out.btm_phi = out.r_btm_zxy_masked.bin.centers;
    out.v_btm_dens_2d = hist3([theta_btm phi_btm],'Nbins',[nbins nbins]);
end
out.data_masked = data_masked;
end