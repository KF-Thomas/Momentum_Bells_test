% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '20191105_halos_3766_shots';
% data_folder = '20191104_halos_attempt_1';
opts.import.dir = fullfile(opts.data_root, data_folder);


opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

opts.num_lim = 1.5e3;

% % Background stuff
cli_header('Setting up for %s', data_folder);
opts.fig_dir = fullfile(this_folder, 'figs', data_folder);
opts.data_src = fullfile(opts.data_root, data_folder);
opts.data_dir = data_folder;
opts.import.cache_save_dir = fullfile(opts.data_root, data_folder, 'cache', 'import\');
opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
opts.index.filename = sprintf('index__%s__%.0f', opts.data_dir);
opts.label = data_folder;
opts.tag = 1;
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

% % import raw data
[data, ~] = import_mcp_tdc_data(opts.import);
%% add labview import
logs = readtable(opts.logfile);
tags = logs{:,5};
%% select a specific shot type if you wish
if opts.tag
shot_type = 'double_halo';
tag_mask = cellfun(@(x) strcmp(x, shot_type), tags');
tag_mask = [tag_mask,zeros(1,length(data.num_counts)-length(tags))];
else
    tag_mask = ones(1,length(data.num_counts));
end
%% set up relevant constants
hebec_constants
%% find centers
opts.cent.visual = 0;
opts.cent.bin_size = 0.5e-5 * [1, 1, 1];
opts.cent.threshold = 2.5;
opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972]}; %time bounds for the different momentum states k=+1,0,-1 respectively
bec = halo_cent(data,opts.cent);
%% run some checks
% atoms number
% laser maybe?
num_check = data.num_counts>opts.num_lim;
is_shot_good = num_check & bec.centre_OK_top' & bec.centre_OK_mid' & bec.centre_OK_btm' & tag_mask;
data_masked = struct_mask(data,is_shot_good);
bec_masked = struct_mask(bec,is_shot_good);
%% convert data to velocity
opts.vel_conv.visual = 1;
opts.vel_conv.title = 'top halo';
% top halo
trap_switch_off= 0;
opts.vel_conv.const.g0 = const.g0;
opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=[0.8,1.2]; %bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = [-0.06,0.06];

opts.vel_conv.bec_center.north = bec_masked.centre_top;
opts.vel_conv.bec_center.south = bec_masked.centre_mid;
opts.vel_conv.bec_width.north = bec_masked.width_top;
opts.vel_conv.bec_width.south = bec_masked.width_mid;
%%
top_halo = halo_vel_conv(data_masked,opts.vel_conv);
%% bottom halo
% if opts.vel_conv.visual
%     stfig('bottom halo')
%     clf
%     xlabel('\(v_x\)')
%     ylabel('\(v_y\)')
%     zlabel('\(v_z\)')
%     hold on
% end
% for this_idx = 1:num_shots % Loop over all shots
%     this_txy = data_masked.counts_txy{this_idx};
%     this_centre = (bec_centre_btm(this_idx, :)+bec_centre_mid(this_idx, :))./2;
%     centred_counts = this_txy - this_centre;
%     txy_btm = bec_centre_btm(this_idx, :)- this_centre;
%     txy_mid = bec_centre_mid(this_idx, :)- this_centre;
% 
%     % Convert to kspace
%     this_outtime = - 0.418707;%this_centre(1)
%     v_btm = txy_to_vel(txy_btm, this_outtime, const.g0, const.fall_distance);
%     v_mid = txy_to_vel(txy_mid, this_outtime, const.g0, const.fall_distance);
%     v_radius = norm(v_btm-v_mid)./2;
%     v_zxy = txy_to_vel(centred_counts, this_outtime, const.g0, const.fall_distance);
%     radius_mask = (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)<v_radius.^2.*v_cut;
%     v_zxy = v_zxy(radius_mask,:);
%     btm_halo.txy{this_idx} = this_txy;
%     btm_halo.vel{this_idx} = v_zxy;
%     clear v_ptr
%     [v_ptr(:,1),v_ptr(:,2),v_ptr(:,3)] = cart2sph(v_zxy(:,2),v_zxy(:,3),v_zxy(:,1));
%     btm_halo.vel_radial{this_idx} = v_ptr;
%     %mask out the BEC's
%     if opts.vel_conv.visual
%         scatter3(v_zxy(:,2),v_zxy(:,3),v_zxy(:,1),'k.')
%     end
% end
%%
v_top_zxy = cell2mat(top_halo.counts_vel');
r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
N_top = top_halo.num_counts;

% v_btm_zxy = cell2mat(btm_halo.vel');
% v_btm = cell2mat(btm_halo.vel_radial');
% M_th = v_btm(:,2) >= theta_range(1) & v_btm(:,2) <= theta_range(2);
% M_ph = (v_btm(:,1) >= min(phi_range) & v_btm(:,1) <= max(phi_range));
% v_btm_masked = v_btm_zxy(M_ph&M_th,:);
% r_dist_btm = v_btm(M_ph&M_th,3);

stfig('radial distribution')
clf
% hist(r_dist_btm,1000)
hold on
hist(r_dist_top,100)
xlabel('r')
ylabel('Freq')
stfig('Counts in halo distribution')
clf
% hist(r_dist_btm,1000)
hold on
hist(N_top,100)
xlabel('N')
ylabel('Freq')
% stfig('top halo masked')
% clf
% scatter3(v_top_masked(:,2),v_top_masked(:,3),v_top_masked(:,1),'k.')
% xlabel('\(v_x\)')
% ylabel('\(v_y\)')
% zlabel('\(v_z\)')
% stfig('bottom halo masked')
% clf
% scatter3(v_btm_masked(:,2),v_btm_masked(:,3),v_btm_masked(:,1),'k.')
% xlabel('\(v_x\)')
% ylabel('\(v_y\)')
% zlabel('\(v_z\)')
%% calculate correlation functions

% back to back (intra halo)

% bb (inter halo)

% cl (intra)
corr_opts.type='1d_cart_cl';%'radial_cl';%'3d_cart_cl';%%
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*6;
one_d_range=0.01;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,50)';
corr_opts.redges=sqrt(linspace(1e-6^2,0.006^2,50));
corr_opts.rad_smoothing=0.0001;

corr_opts.low_mem=nan;
corr_opts.plots=true;
corr_opts.norm_samp_factor=2000;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=true;
corr_opts.sorted_dir=1;
corr_opts.sort_norm=1;
% corr_opts.cl_or_bb = 1;


% corr_opts.one_d_smoothing=nan;
corr_opts.one_d_smoothing=0.00002;
%improved code  
%28.4 s with premask & sort chunks 
%30.12  with premask & no sort chunks 
%28.57 new with no premask
%28.61 s old with no premask

% old code 27.137 s
tic
out=calc_any_g2_type(corr_opts,top_halo.counts_vel);
toc

% cl (inter)