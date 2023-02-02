%% Initializing path
clear all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
%% Import directory
% opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\kieran\Documents\LOCAL-DATA\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';
% opts.data_root = '/Users/tonyyan/Documents/_ANU/_He_BEC_Group/data_root_copy_meh';
opts.data_root = '/Users/tonyyan/Library/CloudStorage/OneDrive-AustralianNationalUniversity/_He_BEC_data_root_copy';


% data_folder = '';%'2023130_new_plates_halo_3_halos';%20221212_new_plates_halo_test_4
data_folder = '20221212_new_plates_halo_test_34_combined_notHe34_just_test_34';
% data_folder = '20221212_new_plates_halo_test_4';

opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;

% opts.import.shot_num = 313:625; %can select specific shots to import

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,4];

opts.num_lim = 1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = -1;%2;%10;%0;% %minimum allowed number in halo 10
opts.halo_N_lim_upper = Inf;%2;%10;%0;% %minimum allowed number in halo 10
y_cut = 11e-3;

z_limits = [-0.9,0.9];%[-0.3,0.3];%[-0.3,0.3];%[-0.4,0.4];%[-0.68,0.68];%[-0.15,0.15];%[-0.15,0.15];%[-0.36,0.36];%
radius_lim = [0.058,0.072];%[0.05,0.07];%[0.03,0.08];%[0.,1.17].*0.065;%[0.79,1.17].*0.065;%[0.61,1.26];%[0.89,1.11];%[0.89,1.16];%[0.9,1.05];%

ang_lim = 45;%35;%angular limit in degrees

plot_dist = true; %do you want to see all the detailed stuff about the halo distributions
opts.corr_center_check = false; %do you want a sceond check

% bec_bounds = {[3.849,3.857],[3.857,3.871],[3.871,3.883],[3.883,3.897]};%{[3.874,3.884],[3.884,3.893]};%{[3.8598,3.871],[3.871,3.8844]};
bec_bounds = {[3.8598,3.871],[3.871,3.8844]};
t0_factor = 3.8772;

do_bb = 1;%if you want to measure the back to back correlations
do_cl = 0;%if you want to measure co linear corrs

%%
%% Background stuff
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];
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
% if ~exist(opts.fig_dir, 'dir')
%     mkdir(opts.fig_dir);
% end
combined_struct = @(S,T) cell2struct(cellfun(@vert_or_horz_cat,struct2cell(S),struct2cell(T),'uni',0),fieldnames(S),1);
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
%data is stored in a structure called data. data.counts_txy{i} is an array
%containing txy data for the ith shot. data.num_counts{i} stores te number
%of counts in the shot. data.shot_num is the shot number. 
% see import_mcp_tdc_data
%% remove any ringing or hotspots
% data_ht_spot=hotspot_mask(data);
% data.counts_txy=data_ht_spot.masked.counts_txy;
% data.num_counts=data_ht_spot.masked.num_counts;
opts.ring_lim = 0.01e-6;%-1;%0.1e-6;%0;%0.101 %how close can points be in time
data_masked = ring_removal(data,opts.ring_lim);

%% add labview import
if opts.tag
    shot_type = 'main';
    if iscell(dirs_list)
        tag_mask = [];
        for ii = 1:length(dirs_list)
            opts.import.dir = dirs_list{ii};
            opts.logfile = fullfile(opts.import.dir, 'log_LabviewMatlab.txt');
            opts_tab = detectImportOptions(opts.logfile);
            opts_tab.Delimiter = {','};
            logs = readtable(opts.logfile,opts_tab);
            tags = logs{:,6};
            %% select a specific shot type if you wish
            tag_mask = [tag_mask,cellfun(@(x) strcmp(x, shot_type), tags')];
            %     tag_mask = [tag_mask,zeros(1,length(data_masked.num_counts)-length(tags))];
        end
        data_masked = struct_mask(data_masked,logical(tag_mask),1);
        
    else
        opts_tab = detectImportOptions(opts.logfile);
        opts_tab.Delimiter = {','};
        logs = readtable(opts.logfile,opts_tab);
        tags = logs{:,6};
        %% select a specific shot type if you wish
        tag_mask = cellfun(@(x) strcmp(x, shot_type), tags');
        tag_mask = [tag_mask,zeros(1,length(data_masked.num_counts)-length(tags))];
        data_masked = struct_mask(data_masked,logical(tag_mask),1);
    end
end

%% get rid of dead shots
num_check = data_masked.num_counts>opts.num_lim;
data_masked = struct_mask(data_masked,num_check);
%% set up relevant constants
hebec_constants

%% find centers
opts.cent.correction = 0;
opts.cent.correction_opts.plots = 0;

% perhaps have these parameters adaptive
opts.cent.visual = 0; %from 0 to 2
opts.cent.savefigs = 0;
opts.cent.threshold = [150,6000,2000.*3].*1e3;  %[130,6000,2000.*3].*1e3;  %[150,80,80].*1e3;  %set in inverse units (Hz for time 1/m for space)
opts.cent.min_threshold = [0,8,8].*1e3;%[16,7,10].*1e3;
opts.cent.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.method = {'margin','average','average'};
% opts.cent.top.method = {'margin','margin','margin'};


opts.cent.t_bounds = bec_bounds; %time bounds for the different momentum states k=-2,-1,0 respectively

for ii = 1:(length(bec_bounds)-1)
    

% opts.cent.t_bounds = {[2.134,2.148],[2.148,2.161],[2.161,2.18],[2.13,2.2]};
% opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% opts.cent.t_bounds = {[1.741,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% t0_factor = 1.7694;

% bec = halo_cent(data_masked,opts.cent);
opts.cent.crop = [opts.cent.t_bounds{ii}; -0.022, 0.011; -0.08, 0.018];
[bec.centre_top,bec.width_top,bec.counts_top,bec.centre_OK_top] =  find_dist_cen(data_masked,opts.cent);

opts.cent.crop = [opts.cent.t_bounds{ii+1}; -0.022, 0.011; -0.08, 0.018];
[bec.centre_btm,bec.width_btm,bec.counts_btm,bec.centre_OK_btm] =  find_dist_cen(data_masked,opts.cent);

%% run some checks
% slosh_cut = ~(abs(bec.centre_mid(:,3))>y_cut);
% num_masked = data_masked.num_counts;
% num_masked(~num_check) = NaN;
% num_outlier = isoutlier(num_masked);
% ~num_outlier &
is_shot_good = bec.centre_OK_top' & bec.centre_OK_btm';%slosh_cut' &  & bec.centre_OK_btm'; & num_check
% data_masked_halo = struct_mask(data_masked,is_shot_good);
bec_masked_halo{ii} = struct_mask(bec,is_shot_good);
data_halo{ii} = struct_mask(data_masked,is_shot_good);
end

%% Find the velocity widths
opts.bec_width.g0 = const.g0;
opts.bec_width.fall_time = 0.417;
% bec_masked_halo = bec_width_txy_to_vel(bec_masked_halo,opts.bec_width);

%% manual centering of halos
% bec_masked_halo.centre_top(:,2:3) = [2.32e-3,-5.2e-3].*ones(size(bec_masked_halo.centre_top(:,2:3)));
% bec_masked_halo.centre_mid(:,2:3) = [2e-3,-5.2e-3].*ones(size(bec_masked_halo.centre_mid(:,2:3)));
% bec_masked_halo.centre_btm = ;

%% convert data to velocity

%% generate top halo
opts.vel_conv.visual = 0;
opts.vel_conv.plot_percentage = 0.95;
opts.vel_conv.const.g0 = const.g0;
opts.vel_conv.const.fall_distance = const.fall_distance;
opts.vel_conv.v_thresh = 0.15; %maximum velocity radius
opts.vel_conv.v_mask=radius_lim;%[0.31,1.56];%[0.89,1.11]; %[0.61,1.26];%bounds on radisu as multiple of radius value
opts.vel_conv.z_mask = z_limits;%[-0.36,0.36];%[-0.65,0.65];%[-0.55,0.55];%[-0.68,0.68]; %[-0.68,0.68]; %in units of radius (standard [-0.76,0.76])
opts.vel_conv.ang_lim = ang_lim; %angular limits of the azimuthal angle
opts.vel_conv.y_mask = [-1.9,1.9];%[-0.8,0.8]; %in units of radius

opts.vel_conv.centering_correction = [0 0 0].*0.5e-3;
opts.vel_conv.phi_correction = [0 0];


% Loop over halos
for ii = 1:(length(bec_bounds)-1)
opts.vel_conv.title = ['halo ', num2str(ii)];

opts.vel_conv.bec_center.north = bec_masked_halo{ii}.centre_top;
opts.vel_conv.bec_center.south = bec_masked_halo{ii}.centre_btm;
opts.vel_conv.bec_width.north = bec_masked_halo{ii}.width_top;
opts.vel_conv.bec_width.south = bec_masked_halo{ii}.width_btm;

% zero velocity point
t0 = ones(size(bec_masked_halo{ii}.centre_top,1),1).*t0_factor;%3.8772;%.*1.7694;%2.1749;%bec_masked_halo.centre_top(:,1);%72;%3.8772
x0 = bec_masked_halo{ii}.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
y0 = bec_masked_halo{ii}.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;

opts.vel_conv.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point

halo_intial = halo_vel_conv(data_halo{ii},opts.vel_conv);

% Mask out halos with nums too low
halo_N_check = halo_intial.num_counts>opts.halo_N_lim & halo_intial.num_counts<opts.halo_N_lim_upper;
halo{ii} = struct_mask(halo_intial,halo_N_check);
bec_halo{ii} = struct_mask(bec_masked_halo{ii},halo_N_check);
end
%% Plot distribution
if plot_dist
for ii = 1:(length(bec_bounds)-1)
v_zxy = cell2mat(halo{ii}.counts_vel_norm);
v_zxy_unnorm = cell2mat(halo{ii}.counts_vel);
r_dist = sqrt(v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2);
r_dist_unnorm = sqrt(v_zxy_unnorm(:,1).^2+v_zxy_unnorm(:,2).^2+v_zxy_unnorm(:,3).^2);
N_current = halo{ii}.num_counts;

stfig('radial distribution');
if ii == 1
clf
else 
    hold on
end
r_hist=smooth_hist(r_dist,'sigma',0.0001);
r_hist_un=smooth_hist(r_dist_unnorm,'sigma',0.0001);
subplot(2,1,1)
hold on
plot(r_hist.bin.centers,r_hist.counts.smooth,'linewidth',1.5)
xlabel('r')
ylabel('Freq')
xlim([min([r_hist.bin.centers]),...
    max([r_hist.bin.centers])])
legend('top','btm')
subplot(2,1,2)
plot(r_hist_un.bin.centers,r_hist_un.counts.smooth,'linewidth',1.5)
hold on
%ylimit = max(r_hist_un.counts.smooth);
%plot([0.130159/2 0.130159/2],[-0.1,ylimit.*2],'k-','linewidth',1.5)
%ylim([0 ylimit.*1.1])
xlabel('r')
ylabel('Freq')
xlim([min(r_hist_un.bin.centers),...
    max(r_hist_un.bin.centers)])
legend('top','btm','expected radius')

stfig('halo comparison');
clf
plot_mask = rand(size(v_zxy,1),1)<0.65;
scatter3(v_zxy(plot_mask,2),v_zxy(plot_mask,3),v_zxy(plot_mask,1),'.')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

stfig(['radial density ',num2str(ii)]);
clf
rad_shift=0.059;

[theta_h,rxy] = cart2pol(v_zxy_unnorm(:,2),v_zxy_unnorm(:,3));
combined_vzr = [v_zxy_unnorm(:,1),rxy;...
        v_zxy_unnorm(:,1),-rxy].*1e3;
    
ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
%caxis([0 4])
    hold on
    plot(zeros(1,1000),linspace(-0.14,0.14,1000).*1e3,'k-','LineWidth',3.8)
    xlabel('$v_r$ (mm/s)')
    %ylabel('$v_z$ (mm/s)')
    colormap('default')
    axis equal
    %caxis([0 9])
    ylim([-0.13,0.12].*1e3)
    ax = gca;
k = 0.02;%
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.

end
end
%% Measure Correlations

%% Measure Squeezing