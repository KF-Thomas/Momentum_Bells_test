% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
data_folder = '20191115_halos_attempt_3';
% data_folder = '20191104_halos_attempt_1';
opts.import.dir = fullfile(opts.data_root, data_folder);


opts.import.force_reimport = true;
opts.import.force_cache_load = ~opts.import.force_reimport;

tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

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

% % import raw data
[data, ~] = import_mcp_tdc_data(opts.import);
%% add labview import
if opts.tag
logs = readtable(opts.logfile);
tags = logs{:,5};
%% select a specific shot type if you wish
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
opts.cent.bin_size = 0.3e-5 * [1, 1, 1];
opts.cent.threshold = 2.5;
opts.cent.t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.8,3.95]}; %time bounds for the different momentum states k=+1,0,-1 respectively
bec = halo_cent(data,opts.cent);
%% run some checks
% atoms number
% laser maybe?
num_check = data.num_counts>opts.num_lim;
is_shot_good = num_check & bec.centre_OK_top' & bec.centre_OK_mid' & bec.centre_OK_btm' & tag_mask;
data_masked = struct_mask(data,is_shot_good);
bec_masked = struct_mask(bec,is_shot_good);
%% Find the velocity widths
opts.bec_width.g0 = const.g0;
opts.bec_width.fall_time = 0.417;
bec_masked = bec_width_txy_to_vel(bec_masked,opts.bec_width);
%% convert data to velocity
opts.vel_conv.visual = 1;
opts.vel_conv.plot_percentage = 0.2;
opts.vel_conv.title = 'top halo';
% top halo
trap_switch_off= 0;
opts.vel_conv.const.g0 = const.g0;
opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=[0.9,1.09]; %bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = [-0.04,0.04];

opts.vel_conv.bec_center.north = bec_masked.centre_top;
opts.vel_conv.bec_center.south = bec_masked.centre_mid;
opts.vel_conv.bec_width.north = bec_masked.width_top;
opts.vel_conv.bec_width.south = bec_masked.width_mid;
%%
top_halo = halo_vel_conv(data_masked,opts.vel_conv);
top_halo.bec_vel_width = (mean(bec_masked.vel_width_top,2)+mean(bec_masked.vel_width_mid,2))./2;% add the average bec width
%% find the mode number
opts.mode_num.qe = 0.08;
top_halo.m = halo_mode_occupancy(top_halo,opts.mode_num);
%% calculated expected correlation amplitude
top_halo.g2 = 1 + 1./top_halo.m;
%% bottom halo

%% plot some histogram checks
v_top_zxy = cell2mat(top_halo.counts_vel');
r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
N_top = top_halo.num_counts;

% v_btm_zxy = cell2mat(btm_halo.vel');
% v_btm = cell2mat(btm_halo.vel_radial');
% M_th = v_btm(:,2) >= theta_range(1) & v_btm(:,2) <= theta_range(2);
% M_ph = (v_btm(:,1) >= min(phi_range) & v_btm(:,1) <= max(phi_range));
% v_btm_masked = v_btm_zxy(M_ph&M_th,:);
% r_dist_btm = v_btm(M_ph&M_th,3);

stfig('radial distribution');
clf
% hist(r_dist_btm,1000)
hold on
hist(r_dist_top,100)
xlabel('r')
ylabel('Freq')
stfig('Counts in halo distribution');
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

            stfig('Halo Num History');
            plot(data_masked.shot_num,...
                top_halo.num_counts,...
                'kx-','LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('N in Halo')
       
            trans_frac = [bec_masked.trans_top';bec_masked.trans_mid';bec_masked.trans_btm';bec_masked.trans_oth'];
            stfig('Transfer Fraction History');
            plot(data_masked.shot_num,...
                trans_frac',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Transfer Fraction')
            legend('top','mid','btm','oth')
            
            stfig('Vel Width History');
            subplot(3,1,1)
            title('Top')
            plot(data_masked.shot_num,...
                bec_masked.vel_width_top','.',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Velocity width (m/s)')
            legend('wz','wx','wy')
            subplot(3,1,2)
            title('Midle')
            plot(data_masked.shot_num,...
                bec_masked.vel_width_mid','.',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Velocity width (m/s)')
            legend('wz','wx','wy')
            subplot(3,1,3)
            title('Bottom')
            plot(data_masked.shot_num,...
                bec_masked.vel_width_btm','.',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Velocity width (m/s)')
            legend('wz','wx','wy')
            
            stfig('Mode Occupancy and Correlation Amplitude')
            subplot(2,1,1)
            plot(data_masked.shot_num,...
                top_halo.m,'kx-',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Mode Occupancy')
            subplot(2,1,2)
            plot(data_masked.shot_num,...
                top_halo.g2,'kx-',...
                'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('$g^2$')

%% calculate correlation functions

% back to back (intra halo)

% bb (inter halo)

% cl (intra)
corr_opts.type='radial_cl';%'1d_cart_cl';%'3d_cart_cl';%%
corr_opts.one_d_dimension=1;
corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*1e-3;
one_d_range=0.06;
corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,150)';
corr_opts.redges=sqrt(linspace(1e-6^2,0.06^2,1500));
corr_opts.rad_smoothing=0.00001;

corr_opts.low_mem=nan;
corr_opts.plots=true;
corr_opts.norm_samp_factor=20;
corr_opts.attenuate_counts=1;
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=0;
corr_opts.sort_norm=0;
% corr_opts.cl_or_bb = 1;


% corr_opts.one_d_smoothing=nan;
corr_opts.one_d_smoothing=0.0001;
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