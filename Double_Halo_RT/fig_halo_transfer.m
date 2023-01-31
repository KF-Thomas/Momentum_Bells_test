%generates figure for halo transfer comparison

%% Initializing path
% clear all;
 close all;
this_folder = fileparts(which(mfilename));
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')
%% Import directory
opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
% opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';

%Z:\EXPERIMENT-DATA\2020_Momentum_Bells\pulse_characterisation\20210129_double_mirror
%Z:\EXPERIMENT-DATA\2020_Momentum_Bells\pulse_characterisation\20210201_bragg_pulse_testing\double_mirror

% data_folder = 'single_halo_data\20210316_k=-1,-2_various_tests\20210316_k=-1,-2_mirror';
% data_folder_norm = 'single_halo_data\20210316_k=-1,-2_various_tests\20210316_k=-1,-2_test_1';

data_folder='20221102_new_plates_halo_test';%'20211206_scaning_across_freq\bs_high_100_kHz_evap_0_844_MHz';%'20211206_scaning_across_freq\double_mirror_120_kHz_and_80_kHz';%
% data_folder='20211206_scaning_across_freq\mr_98_kHZ';
% data_folder='20211206_scaning_across_freq\changing_width_and_amp_2';
data_folder_norm='20221102_new_plates_halo_test';%'20211206_scaning_across_freq\norm_0_844_MHz';

% data_folder = '20220306_k=-1,-2_halo_mirror';
% data_folder_norm = '20220306_k=-1,-2_halo_norm';

% data_folder='20220210_k=-1,-2_pulse_characteristic\k=-1,-2_halo_mirror';%'20211206_scaning_across_freq\double_mirror_120_kHz_and_80_kHz';%
% 
% data_folder_norm='20220210_k=-1,-2_pulse_characteristic\k=-1,-2_halo';

%double halo
% data_folder='20211027_validating_bragg_pulses\source pulse';
% data_folder='20211206_scaning_across_freq\norm_0_839_MHz';

%single halo (top)
% data_folder='20211206_scaning_across_freq\norm_0_844_MHz';

opts.import.dir = fullfile(opts.data_root, data_folder);
opts.import.force_reimport = false;
opts.import.force_cache_load = ~opts.import.force_reimport;


%% Chose which halo(s) to analyse
opts.do_top_halo = 1;% analyse the top halo?
opts.do_btm_halo = 1;% analyse the bottom halo?

%% Chose if you want to look at a narrow or wide slice of the halo
slice_type = 'medium';
if strcmp(slice_type,'narrow')
    opts.vel_conv.top.z_mask = [-0.4,0.4];%
    opts.vel_conv.btm.z_mask = [-0.4,0.4];%in units of radius ([-0.68,0.68])
elseif strcmp(slice_type,'extra wide')
    opts.vel_conv.top.z_mask = [-0.87,0.87];%
    opts.vel_conv.btm.z_mask = [-0.87,0.87];%in units of radius ([-0.68,0.68])
elseif strcmp(slice_type,'medium')
    opts.vel_conv.top.z_mask = [-0.6,0.6];%
    opts.vel_conv.btm.z_mask = [-0.6,0.6];%in units of radius ([-0.68,0.68])
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
% 
opts.vel_conv.top.z_mask = [-5,3];%[-3,1];%[-0.9,0.9];%[-1,1];%[-0.9,0.9];
opts.vel_conv.btm.z_mask = [-0.9,0.9];%[-0.01,0.02];%[-1,1];%[-0.9,0.9];%in units of radius ([-0.68,0.68])

opts.radius_lim = [0,3.5];%[-Inf,Inf];%[-0.03,0.08];%[0.00,0.14];%[0.01.*0.065,0.08];%[0.01.*0.065,0.1];%[0.05,0.07];%[0.3,1.61].*0.065;%[0.61,1.26].*0.065;%[0.9,1.05];%[0.89,1.11];
opts.ang_lim = 120;%92;%angular limit in degrees

%% Import parameters
tmp_xlim=[-35e-3, 35e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-35e-3, 35e-3];
tlim=[0,6];
opts.import.txylim=[tlim;tmp_xlim;tmp_ylim];

opts.num_lim = 0e3;%2.1e3;%0.5e3;% %minimum atom number 1.5e3
opts.halo_N_lim = 0;%2;%10;%0;% %minimum allowed number in halo 10

opts.plot_dist = false; %do you want to see all the detailed stuff about the halo distributions

opts.cent.nan_cull = false; %do you want to cull nan's

%% find centers
opts.cent.visual = 0; %from 0 to 2
opts.cent.savefigs = 0;
opts.cent.correction = 0;%2
opts.cent.correction_opts.plots = 0;

opts.cent.top.visual = 0; %from 0 to 2
opts.cent.top.savefigs = 0;
opts.cent.top.threshold = [130,5000,5000].*1e3;%130
opts.cent.top.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[16,7,10].*1e3;
opts.cent.top.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.top.method = {'margin','average','average'};

opts.cent.mid.visual = 0; %from 0 to 2
opts.cent.mid.savefigs = 0;
opts.cent.mid.threshold = [130,5000,5000].*1e3;
opts.cent.mid.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[16,7,10].*1e3;
opts.cent.mid.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.mid.method = {'margin','average','average'};

opts.cent.btm.visual = 0; %from 0 to 2
opts.cent.btm.savefigs = 0;
opts.cent.btm.threshold = [130,5000,5000].*1e3;%[130,2000,2000].*1e3;
opts.cent.btm.min_threshold = [10,3,3].*1e3;%[16,3,3].*1e3;%[0,0,0].*1e3;%[16,13,13].*1e3;%[16,7,10].*1e3;
opts.cent.btm.sigma = [6.7e-5,16e-5,16e-5];%[8e-5,25e-5,25e-5];
opts.cent.btm.method = {'margin','average','average'};

% opts.cent.t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% opts.cent.t_bounds = {[1.741,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};
% opts.cent.t_bounds = {[2.134,2.148],[2.148,2.161],[2.161,2.18],[2.13,2.2]};
 opts.cent.t_bounds = {[3.844,3.8598],[3.8598,3.871],[3.871,3.8844],[3.75,4]};%time bounds for the different momentum states
% opts.cent.t_bounds = {[5.350,5.356],[5.361,5.367],[5.372,5.380],[5.34,5.39]};%time bounds for the different momentum states (for full evap settings)
%%
out = halo_func(data_folder,opts);
out_norm = halo_func(data_folder_norm,opts);

if opts.do_top_halo && opts.do_btm_halo
v_top_dens = out.v_top_dens;
cal_dens_top = 1.*out_norm.v_top_dens;
v_btm_dens = out.v_btm_dens;
cal_dens_btm = 1.*out_norm.v_btm_dens;
trans_ratio_top = (v_top_dens(:,2)-cal_dens_top(:,2))./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_top(:,2).*2);
trans_ratio_btm = (v_btm_dens(:,2)-cal_dens_btm(:,2))./(v_top_dens(:,2)+v_btm_dens(:,2)-cal_dens_btm(:,2).*2);
end
%%
transfer_func = @(a,b,c) (a-c)./(a+b-c.*2);
cal_dens = 1.*out_norm.v_dens;
out_dens = 1.*out.v_dens;
bin_num = length(cal_dens);
rad_val = 6.50795e-02;%
[rad_tol,rad_indx]=min(abs(rad_val.*2-(out.v_dens_cen(2:end)-out.v_dens_cen(1))));
cal_dens_b = cal_dens(1:bin_num-rad_indx);
cal_dens_t = cal_dens(1+rad_indx:end);
dens_b = out_dens(1:bin_num-rad_indx);
dens_t = out_dens(1+rad_indx:end);
norm_C =  cal_dens_b+cal_dens_t./(dens_b+dens_t);
% dens_b = norm_C.*dens_b;
% dens_t = norm_C.*dens_t;


c1 = transfer_func(dens_t,dens_b,cal_dens_t);
cm1 = transfer_func(dens_b,dens_t,cal_dens_b);
[a,b,c,x] = transfer_perc(out_dens,cal_dens,rad_indx,out.v_dens_cen);
figure
subplot(2,1,1)
% plot(x,b)
% hold on
% plot(x,c)
plot(x,a)
hold on
plot(x,c)
ylim([0 1])
subplot(2,1,2)
plot(out.v_dens_cen(:),out_dens(:))
hold on
plot(out.v_dens_cen(:),cal_dens(:))

%% generate transfer figure
stfig('transfer coefficents')
clf
plot(sin(out.top_phi)+1,trans_ratio_top,'r','LineWidth',1.5)
hold on
ylim([0 1])
plot(sin(out.btm_phi)-1,trans_ratio_btm,'b--','LineWidth',1.5)
plot(sin(out.top_phi)+1,1-trans_ratio_top,'k-.','LineWidth',1.5)
plot(sin(out.btm_phi)-1,1-trans_ratio_btm,'k-.','LineWidth',1.5)



%% full 3d comparison
plt_p = 1;
rad_shift=0.059;
z_shift_top = [rad_shift.*ones(size(out.v_top_zxy_unnorm,1),1),zeros(size(out.v_top_zxy_unnorm,1),2)];
z_shift_btm = [rad_shift.*ones(size(out.v_btm_zxy_unnorm,1),1),zeros(size(out.v_btm_zxy_unnorm,1),2)];

z_shift_top_norm = [rad_shift.*ones(size(out_norm.v_top_zxy_unnorm,1),1),zeros(size(out_norm.v_top_zxy_unnorm,1),2)];
z_shift_btm_norm = [rad_shift.*ones(size(out_norm.v_btm_zxy_unnorm,1),1),zeros(size(out_norm.v_btm_zxy_unnorm,1),2)];
    

%%
    stfig('density of halos radius');
    clf
    set(gcf,'position',[516    42  1600   774])%[916    42   389   774]
    h(3)=subplot(1,3,3);
plot(a,(x+rad_shift).*1e3,'b','LineWidth',3)
hold on
%plot(1-a-temp,(x+rad_shift).*1e3,'k--','LineWidth',1.5)
% plot(temp,(x+rad_shift).*1e3,'r','LineWidth',1.5)
ylim([-0.13,0.12].*1e3)
xlim([0 1])
set(gca,'FontSize',25)
set(gca,'yticklabel',[])
set(gca, 'Layer','top')
xlabel('Proportion')
legend('$|C_{-1}|^2$','$|C_{+1}|^2$')
ax = gca;
properties(ax)
k = 0.02;%
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.

% h(1)=subplot(2,2,1);
all_counts=cell2mat(out.data_masked.counts_txy.');
% % all_counts=all_counts((all_counts(:,1)<1.8)&(all_counts(:,1)>1.7),:);
% smooth_counts=smooth_hist(all_counts(:,1),'sigma',1e-5,'lims',[1.74,1.78],'bin_num',10000);
% plot(smooth_counts.bin.centers,smooth_counts.count_rate.smooth./109./1e3,'k-')
% ylabel('Flux (kHz/shot)')
% xlabel('Arrival time (s)')
% set(h(1), 'position', [0.3, 0.7, 0.5, 0.2] );
% set(h(1),'FontSize',13)
h(2)=subplot(1,3,2);
if opts.do_top_halo && opts.do_btm_halo
%     stfig('density of halos radius');
%     clf
%     combined_vzr = [out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),out.rxy_top;...
%         out.v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),out.rxy_btm
%         out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),-out.rxy_top;...
%         out.v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),-out.rxy_btm].*1e3;
combined_vzr = [out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),out.rxy_top;...
        out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),-out.rxy_top].*1e3;
    
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
    set(gca,'yticklabel',[])
    set(gca, 'Layer','top')
    ax = gca;
properties(ax)
k = 0.02;%
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.
    
    
elseif opts.do_top_halo

    
    combined_vzr = [v_top_zxy_unnorm(:,1)+z_shift_top(:,1),rxy_top;...
        v_top_zxy_unnorm(:,1)+z_shift_top(:,1),-rxy_top].*1e3;
    ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
    hold on
    plot(zeros(1,1000),linspace(-0.14,0.14,1000).*1e3,'k-','LineWidth',3.8)
    xlabel('$v_r$ (mm/s)')
    ylabel('$v_z$ (mm/s)')
    colormap('default')
    axis equal
    %caxis([0 9])
    ylim([-0.13,0.12].*1e3)
elseif opts.do_btm_halo
%     stfig('density of halos radius');
%     clf
    combined_vzr = [v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),rxy_btm];
    ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
    xlabel('$v_r$ (m/s)')
    ylabel('$v_z$ (m/s)')
    colormap('default')
    axis equal
    %caxis([0 9])
    
end
caxis([0 12])

box on
xlim([-73 73])
set(gca,'FontSize',25)
c=colorbar('northoutside');
c.Label.String = 'Counts density (counts/(mm/s)$^2$/shot)';
c.Label.Interpreter = 'latex';
c.FontName = 'Times';

h(1)=subplot(1,3,1);

%     stfig('density of halos radius');
%     clf
%     combined_vzr = [out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),out.rxy_top;...
%         out.v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),out.rxy_btm
%         out.v_top_zxy_unnorm(:,1)+z_shift_top(:,1),-out.rxy_top;...
%         out.v_btm_zxy_unnorm(:,1)-z_shift_btm(:,1),-out.rxy_btm].*1e3;
combined_vzr = [out_norm.v_top_zxy_unnorm(:,1)+z_shift_top_norm(:,1),out_norm.rxy_top;...
        out_norm.v_top_zxy_unnorm(:,1)+z_shift_top_norm(:,1),-out_norm.rxy_top].*1e3;
    
[o1,o2,o3]=ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
    hold on
    plot(zeros(1,1000),linspace(-0.14,0.14,1000).*1e3,'k-','LineWidth',3.8)
    xlabel('$v_r$ (mm/s)')
    ylabel('$v_z$ (mm/s)')
    colormap('default')
    axis equal

    ylim([-0.13,0.12].*1e3)

caxis([0 12])
box on
xlim([-73 73])
set(gca,'FontSize',25)
set(gca, 'Layer','top')
ax = gca;
properties(ax)
k = 0.02;%
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.



set(h(3), 'position', [0.67, 0.15, 0.2, 0.7] );
set(h(2), 'position', [0.27, 0.15, 0.48, 0.7] );
set(h(1), 'position', [0.00, 0.15, 0.48, 0.7] );     



% Get Position of the figure 1
%pos = get(gcf, 'Position'); % gives x left, y bottom, width, height
%x = pos(1); y = pos(2); w = pos(3); h = pos(4);

%%
function [cm1,c0,c1,x]=transfer_perc(c,cal,rad_indx,bin_cen)
transfer_func = @(a,b,c) (a-c)./(a+b-c.*2);
bin_num = length(c);
% cal = cal;
m_index = 1:(bin_num-2*rad_indx);
z_index = rad_indx+1:bin_num-rad_indx;
p_index = 2*rad_indx+1:bin_num;
C_norm = (cal(m_index)+cal(z_index))./(c(m_index)+c(z_index));
cm1 = transfer_func(C_norm.*c(m_index),C_norm.*c(z_index),cal(m_index));
cm1(cm1<0)=0;
cm1(cm1>1)=0;
c1 = transfer_func(c(p_index),c(z_index),cal(p_index));
c1(c1<0)=0;
c1(c1>1)=0;
c_sum = cm1+c1;
c1_comp= (c(p_index)-cal(p_index))./(cal(z_index)-cal(p_index));

c0 = 1-(cm1+c1);
N_val = (cm1+c0+c1);
c1 = c1./N_val;
c0 = c0./N_val;
cm1 = cm1./N_val;

x = bin_cen(z_index);

end
