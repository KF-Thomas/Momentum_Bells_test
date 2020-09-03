%% Initializing path
% clear all;
% close all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
%% Import directory
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '';
% data_folder = '20200807_k=0,-1,-2_halos_data_2';
% data_folder = '20200824_k=0,-1_splitter_attempt';
% data_folder = '20200824_k=0,-1_pre_mirror';
% data_folder = '20200825_attempting_different_mirror_settings';
% data_folder = '20200824_k=0,-1_mirror_attempt_1';
% data_folder = '20200825_k=0,-1_mirror_attempt_2';
% data_folder = '20200828_k=0,-1_mirror_attempt_3';
% data_folder = '20200831_k=0,-1_mirror_attempt_4';
% data_folder = '20200831_k=0,-1_mirror_attempt_5';
% data_folder = '20200831_k=0,-1_mirror_attempt_6';
% data_folder = '20200831_k=0,-1_mirror_attempt_7';
% data_folder = '20200901_k=0,-1_mirror_attempt_8';
% data_folder = '20200901_k=0,-1_mirror_vs_amp\Pamp_0_5';
data_folder = '20200901_k=0,-1_mirror_vs_amp\Pamp_10';
data_folder = '20200901_k=0,-1_mirror_vs_amp\Pamp_0_25';
data_folder = '20200901_k=0,-1_mirror_vs_amp\Pamp_0_25';
opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

%% Chose which halo(s) to analyse
opts.do_top_halo = 1;% analyse the top halo?
opts.do_btm_halo = 1;% analyse the bottom halo?

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 2.5e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10

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
opts.cent.top.threshold = [130,2000,2000].*1e3;
opts.cent.top.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.top.method = {'margin','average','average'};

opts.cent.mid.visual = 0; %from 0 to 2
opts.cent.mid.savefigs = 0;
opts.cent.mid.threshold = [130,2000,2000].*1e3;
opts.cent.mid.min_threshold = [16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.mid.method = {'margin','average','average'};

opts.cent.btm.visual = 2; %from 0 to 2
opts.cent.btm.savefigs = 0;
opts.cent.btm.threshold = [15,5000,5000].*1e3;%[130,2000,2000].*1e3;
opts.cent.btm.min_threshold = [16,7,10].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
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
opts.vel_conv.top.z_mask = [-0.82,0.82]; %in units of radius (standard [-0.76,0.76] or [-0.68,0.68])
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
opts.vel_conv.btm.z_mask = [-0.82,0.82]; %in units of radius ([-0.68,0.68])
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

%     top_halo = halos.top_halo;
%     bottom_halo = halos.bottom_halo;
%     bec_masked = halos.bec;

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

stfig('average radius vs angle')
clf
nbins=50;
theta_bins = linspace(-pi,pi,nbins);
phi_bins = linspace(-pi/2,pi/2,nbins);
if opts.do_top_halo
[theta_top,~] = cart2pol(v_top_zxy(:,2),v_top_zxy(:,3));
phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
else
    theta_top = [];
    phi_top = [];
end

if opts.do_btm_halo
[theta_btm,~] = cart2pol(v_btm_zxy(:,2),v_btm_zxy(:,3));
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
    v_btm_dens(ii,1) = size(r_btm_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
    v_top_dens(ii,1) = size(r_top_zxy_masked,1)/(theta_bins(ii+1)-theta_bins(ii));
    
    r_btm_zxy_masked = r_dist_btm_unnorm(phi_bins(ii)<phi_btm & phi_btm<=phi_bins(ii+1));
    r_top_zxy_masked = r_dist_top_unnorm(phi_bins(ii)<phi_top & phi_top<=phi_bins(ii+1));
    v_btm_r(ii,2) = mean(r_btm_zxy_masked);
    v_top_r(ii,2) = mean(r_top_zxy_masked);
    v_btm_dens(ii,2) = size(r_btm_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));
    v_top_dens(ii,2) = size(r_top_zxy_masked,1)/(phi_bins(ii+1)-phi_bins(ii));
    
    theta(ii) = mean(theta_bins(ii:(ii+1)));
    phi(ii) = mean(phi_bins(ii:(ii+1)));
end
v_btm_dens = v_btm_dens./size(bottom_halo.counts_vel,1);
v_top_dens = v_top_dens./size(top_halo.counts_vel,1);

subplot(2,1,1)
plot(theta,v_btm_r(:,1),'linewidth',1.5)
hold on
plot(theta,v_top_r(:,1),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\theta\)')
subplot(2,1,2)
plot(phi,v_btm_r(:,2),'linewidth',1.5)
hold on
plot(phi,v_top_r(:,2),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\phi\)')

stfig('density vs angle')
clf
subplot(2,1,1)
plot(theta,v_btm_dens(:,1),'linewidth',1.5)
hold on
plot(theta,v_top_dens(:,1),'linewidth',1.5)
legend('bottom','top')
ylabel('Average density')
xlabel('\(\theta\)')
subplot(2,1,2)
plot(phi,v_btm_dens(:,2),'linewidth',1.5)
hold on
plot(phi,v_top_dens(:,2),'linewidth',1.5)
legend('bottom','top')
ylabel('Average density')
xlabel('\(\phi\)')
if opts.do_btm_halo && opts.do_top_halo
stfig('density ratio vs angle');
clf
plot(phi,(v_btm_dens(:,2))./(v_top_dens(:,2)),'linewidth',1.5)
ylabel('density ratio')
xlabel('\(\phi\)')

stfig('density visibility vs angle');
clf
plot(phi,(v_btm_dens(:,2)-v_top_dens(:,2))./(v_btm_dens(:,2)+v_top_dens(:,2)),'linewidth',1.5)
ylabel('density vis')
xlabel('\(\phi\)')

end

stfig('spherical density plot')
for ii = 1:(nbins-1)
    for jj = 1:(nbins-1)
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
