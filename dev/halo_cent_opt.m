function out = halo_cent_opt(data_masked_halo,bec_masked_halo,correction,top_or_btm)
hebec_constants
opts.halo_N_lim = 0;
%% BACK TO BACK (in the same halo)
%% general settings
corr_opts.verbose = false;
corr_opts.print_update = false;
corr_opts.timer=false;

corr_opts.plots = false;
corr_opts.fit = false;
corr_opts.calc_err = false;

% variables for calculating the error
corr_opts.samp_frac_lims=[0.65,0.9];
corr_opts.num_samp_frac=5;
corr_opts.num_samp_rep=5;

corr_opts.attenuate_counts=1;

corr_opts.type='radial_bb';%'1d_cart_bb';%
% corr_opts.one_d_dimension = 2;
% corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*3e-3;
one_d_range=0.017;%0.02
% one_d_range=0.16;%0.02
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,75));%100 or 80 or 85 or 95
% corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,150);
corr_opts.rad_smoothing=nan;
corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

corr_opts.norm_samp_factor=1500;%1500;
corr_opts.sample_proportion=1.0;%0.65;%1500;
corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=nan;
corr_opts.sort_norm=0;

corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit


%% convert data to velocity
% zero velocity point
t0 = ones(size(bec_masked_halo.centre_top,1),1).*3.8772;%bec_masked_halo.centre_top(:,1);%72;%
x0 = bec_masked_halo.centre_top(:,2);%ones(size(bec_masked_halo.centre_top,1),1).*-0.0041;%%-0.00444892593829574;
y0 = bec_masked_halo.centre_top(:,3);%ones(size(bec_masked_halo.centre_top,1),1).*0.0078;%0.00645675151404596;

switch top_or_btm
    case 1        
        %% generate top halo
        opts.vel_conv.top.visual = 0;
        opts.vel_conv.top.plot_percentage = 0.95;
        opts.vel_conv.top.title = 'top halo';
        opts.vel_conv.top.const.g0 = const.g0;
        opts.vel_conv.top.const.fall_distance = const.fall_distance;
        opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
        opts.vel_conv.top.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
        opts.vel_conv.top.z_mask = [-0.65,0.65];%[-0.68,0.68]; %[-0.68,0.68]; %in units of radius (standard [-0.76,0.76])
        opts.vel_conv.top.y_mask = [-0.8,0.8]; %in units of radius
        opts.vel_conv.top.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point
        
        opts.vel_conv.top.centering_correction = correction.*0.5e-3;%[-0.2223,0.5662,-0.8083].*0.5e-3;%[0,0,0]; %[-0.73,0.822,-1.209].*0.5e-3;%[-0.96226,0.788847,-1.11].*0.5e-3;%[-0.6519,0.7836,-1.167].*0.5e-3;%[0,0,0]; %[0.677,0.9842,-1.139].*0.5e-3;%[0,0,0]; %;[3.145e-1,1.313,-1.1705].*0.5e-3;%[0,0,0]; %[0,0,0]; %correctoin shift to the centering in m/s
        opts.vel_conv.top.phi_correction = [0 0];
        
        opts.vel_conv.top.bec_center.north = bec_masked_halo.centre_top;
        opts.vel_conv.top.bec_center.south = bec_masked_halo.centre_mid;
        opts.vel_conv.top.bec_width.north = bec_masked_halo.width_top;
        opts.vel_conv.top.bec_width.south = bec_masked_halo.width_mid;
        
        %%
        top_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.top);
        
        %% mask out halos with nums to low
        halo_N_check = top_halo_intial.num_counts>opts.halo_N_lim;
        top_halo = struct_mask(top_halo_intial,halo_N_check);
        
        %% TOP HALO BACK TO BACK
        
        corr_opts.fig = 'top halo bb corr';
        corrs = calc_any_g2_type(corr_opts,top_halo.counts_vel');
        %     corrs=calc_any_g2_type(corr_opts,top_halo.counts_vel_norm');
        out = corrs.norm_g2.g2_amp(1);%.g2peak;
        
        
    case 0
        %% generate bottom halo
        opts.vel_conv.btm.visual = 0;
        opts.vel_conv.btm.plot_percentage = 0.95;
        opts.vel_conv.btm.title = 'bottom halo';
        opts.vel_conv.btm.const.g0 = const.g0;
        opts.vel_conv.btm.const.fall_distance = const.fall_distance;
        opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
        opts.vel_conv.btm.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
        opts.vel_conv.btm.z_mask = [-0.65,0.65];%[-0.68,0.68]; %[-0.68,0.68]; %in units of radius
        opts.vel_conv.btm.y_mask = [-0.8,0.8]; %in units of radius
        opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point
        
        opts.vel_conv.btm.centering_correction = correction.*0.5e-3;%[-0.1733,1.075,-0.9288].*0.5e-3; %[0,0,0].*0.5e-3;%[-0.1733,1.075,-0.9288].*0.5e-3; %[0.1169,1.606,-1.438].*0.5e-3;%[[0.205,1.7893,-1.4207].*0.5e-3;%[0,0,0]; %[-0.452,1.76,-1.561].*0.5e-3;%[-0.1762,1.6035,-1.029].*0.5e-3;%[0,1.73,-1.45].*0.5e-3; %[-2,-2,1.5].*-0.5e-3; %correctoin shift to the centering in m/s
        opts.vel_conv.btm.phi_correction = [0 0];
        
        opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
        opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
        opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
        opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;
        
        %%
        bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);
        
        %% mask out halos with nums to low
        halo_N_check = bottom_halo_intial.num_counts>opts.halo_N_lim;
        bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
        
        %% BOTTOM HALO BACK TO BACK
        
        corr_opts.fig='bottom halo bb corr';
        corrs=calc_any_g2_type(corr_opts,bottom_halo.counts_vel');
        %     corrs=calc_any_g2_type(corr_opts,bottom_halo.counts_vel_norm');
        out = corrs.norm_g2.g2_amp(1);%.g2peak;
    case 2
        %% generate top halo
        opts.vel_conv.top.visual = 0;
        opts.vel_conv.top.plot_percentage = 0.95;
        opts.vel_conv.top.title = 'top halo';
        opts.vel_conv.top.const.g0 = const.g0;
        opts.vel_conv.top.const.fall_distance = const.fall_distance;
        opts.vel_conv.top.v_thresh = 0.15; %maximum velocity radius
        opts.vel_conv.top.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
        opts.vel_conv.top.z_mask = [-0.55,0.55];%[-0.68,0.68]; %[-0.68,0.68]; %in units of radius (standard [-0.76,0.76])
        opts.vel_conv.top.y_mask = [-1.9,1.9]; %in units of radius
        opts.vel_conv.top.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%%bec_masked_halo.centre_top;%bec_masked_halo.centre_mid; %use the mid BEC as the zero momentum point
        
        opts.vel_conv.top.centering_correction = [-0.40102     0.39517     0.40242].*0.5e-3;%correction(1:3).*0.5e-3;%[-0.2223,0.5662,-0.8083].*0.5e-3;%[0,0,0]; %[-0.73,0.822,-1.209].*0.5e-3;%[-0.96226,0.788847,-1.11].*0.5e-3;%[-0.6519,0.7836,-1.167].*0.5e-3;%[0,0,0]; %[0.677,0.9842,-1.139].*0.5e-3;%[0,0,0]; %;[3.145e-1,1.313,-1.1705].*0.5e-3;%[0,0,0]; %[0,0,0]; %correctoin shift to the centering in m/s
        opts.vel_conv.top.phi_correction = [0 0];
        
        opts.vel_conv.top.bec_center.north = bec_masked_halo.centre_top;
        opts.vel_conv.top.bec_center.south = bec_masked_halo.centre_mid;
        opts.vel_conv.top.bec_width.north = bec_masked_halo.width_top;
        opts.vel_conv.top.bec_width.south = bec_masked_halo.width_mid;
        
        %%
        top_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.top);
        
        %% mask out halos with nums to low
        halo_N_check = top_halo_intial.num_counts>opts.halo_N_lim;
        top_halo = struct_mask(top_halo_intial,halo_N_check);
        
     
        %% generate bottom halo
        opts.vel_conv.btm.visual = 0;
        opts.vel_conv.btm.plot_percentage = 0.95;
        opts.vel_conv.btm.title = 'bottom halo';
        opts.vel_conv.btm.const.g0 = const.g0;
        opts.vel_conv.btm.const.fall_distance = const.fall_distance;
        opts.vel_conv.btm.v_thresh = 0.15; %maximum velocity radius
        opts.vel_conv.btm.v_mask=[0.89,1.11]; %bounds on radisu as multiple of radius value
        opts.vel_conv.btm.z_mask = [-0.55,0.55];%[-0.68,0.68]; %[-0.68,0.68]; %in units of radius
        opts.vel_conv.btm.y_mask = [-1.9,1.9]; %in units of radius
        opts.vel_conv.btm.center = [t0,x0,y0];%bec_masked_halo.centre_top;%ones(size(bec_masked_halo.centre_top,1),1).*[t0,x0,y0];%,bec_masked_halo.centre_top; %use the mid BEC as the zero momentum point
        
        opts.vel_conv.btm.centering_correction = [2.3195      1.2916     -1.5342].*0.5e-3;%correction(4:6).*0.5e-3;%[-0.1733,1.075,-0.9288].*0.5e-3; %[0,0,0].*0.5e-3;%[-0.1733,1.075,-0.9288].*0.5e-3; %[0.1169,1.606,-1.438].*0.5e-3;%[[0.205,1.7893,-1.4207].*0.5e-3;%[0,0,0]; %[-0.452,1.76,-1.561].*0.5e-3;%[-0.1762,1.6035,-1.029].*0.5e-3;%[0,1.73,-1.45].*0.5e-3; %[-2,-2,1.5].*-0.5e-3; %correctoin shift to the centering in m/s
        opts.vel_conv.btm.phi_correction = correction(1:2);
        
        opts.vel_conv.btm.bec_center.north = bec_masked_halo.centre_mid;
        opts.vel_conv.btm.bec_center.south = bec_masked_halo.centre_btm;
        opts.vel_conv.btm.bec_width.north = bec_masked_halo.width_mid;
        opts.vel_conv.btm.bec_width.south = bec_masked_halo.width_btm;
        
        %%
        bottom_halo_intial = halo_vel_conv(data_masked_halo,opts.vel_conv.btm);
        
        %% mask out halos with nums to low
        halo_N_check = bottom_halo_intial.num_counts>opts.halo_N_lim;
        bottom_halo = struct_mask(bottom_halo_intial,halo_N_check);
        
        %% BETWEEN HALO BACK TO BACK
        corr_opts.fig='between halo bb corr';
        corr_opts.type='radial_bb';%'3d_cart_bb';%'1d_cart_bb';%
        corr_opts.one_d_dimension=3;
        corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*0.05;
        % one_d_range=0.16;
        % corr_opts.one_d_window=[[-1,1];[-1,1];[-1,1]]*5e-2;
        one_d_range=0.02;
        corr_opts.one_d_edges=linspace(-one_d_range,one_d_range,100)';
        corr_opts.redges=sqrt(linspace(0,one_d_range^2,40));%(5e-3)^2
        corr_opts.rad_smoothing=nan;
        corr_opts.direction_labels = {'z','x','y'};
        corr_opts.low_mem=true;
        corr_opts.sampling_method='complete';%
        corr_opts.norm_samp_factor=1500;%1500;
%         corr_opts.sample_proportion=global_sample_portion;%1500;
        corr_opts.norm_samp_factor=1;
        corr_opts.attenuate_counts=1;
        corr_opts.do_pre_mask=false;
        corr_opts.sorted_dir=1;
        corr_opts.sort_norm=1;
        corr_opts.progress_updates=5;
        
        corr_opts.gaussian_fit = false; %ensure it always uses a gaussian fit
        
        % both_halo_counts = [top_halo.counts_vel_norm';bottom_halo.counts_vel_norm'];
        both_halo_counts = [top_halo.counts_vel';bottom_halo.counts_vel'];
        
        corr_opts.one_d_smoothing=nan;
        % corr_opts.one_d_smoothing=0.0008;
        
        corrs=calc_any_g2_type(corr_opts,both_halo_counts);
        out = corrs.norm_g2.g2_amp(1);%.g2peak;
end