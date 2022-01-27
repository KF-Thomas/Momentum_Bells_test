%data vs dk
%% 
%% Setting up variables
%     halos.top_halo = top_halo;
%     halos.bottom_halo = bottom_halo;
%     halos.bec = bec_halo;

for ii = 1:length(phi_vec)      
    %% Separate data into the four ports
    ports = {};
    [ports.top_left, ports.top_right] = separate_ports(out_data{ii}.top_halo,0);
    [ports.bottom_left, ports.bottom_right] = separate_ports(out_data{ii}.bottom_halo,0);
    
    %% Quantum correlator E
    opts_E.calc_err = true;
    opts_E.plots = false;
    opts_E.verbose = false;
    opts_E.fit = false;
    opts_E.norm = false; %use normalised or unnormalised data
%     opts_E.sample_proportion = 1.0;%0.1;
    lambda = 0.1;
%     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
    opts_E.delta_kd = [1e-3,1.5e-3.*4,6e-3.*4];%
    opts_E.dim = 1;
    opts_E.sample_proportion = 1.0;
    opts_E.num_samp_rep = 50;
    
    %% BACK TO BACK (in the same halo)
    corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    
    global_sample_portion = 1.0;
    
    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    dkr = (dkx.*dky.*dkz).^(1/3);
    
    % variables for calculating the error
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.type='1d_cart_bb';%'radial_bb';%
    corr_opts.plots = true;
    corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
    corr_opts.fit = false;
    corr_opts.calc_err = do_g2_err;
    corr_opts.one_d_dimension = 2;
    corr_opts.two_d_dimensions = [2,3];
    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
    
    one_d_range=0.05;

    num_pts_cart = round(one_d_range./opts_E.delta_kd(corr_opts.one_d_dimension));
    num_pts_rad = round(one_d_range./dkr);
    num_pts_1 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(1)));
    num_pts_2 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(2)));
    corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
    corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);
    corr_opts.two_d_edges = {linspace(-one_d_range,one_d_range,num_pts_1)',linspace(-one_d_range,one_d_range,num_pts_2)'};
    corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';
    
    corr_opts.rad_smoothing=nan;
    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;
    
    corr_opts.norm_samp_factor=1500;%1500;
    corr_opts.sample_proportion=1.0;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    
    corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
    
    %% TOP HALO BACK TO BACK

        
        [E_val(ii), corrs.ports] = E(ports,opts_E);
        
        out_corrs{ii} = corrs.ports;
end

num_phi = length(out_corrs);

for ii = 1:num_phi
    [k_z_12{ii},k_z{ii},w_z_12{ii}] = g2_cumulative_avg(out_corrs{ii}.g12,0.1);
    [k_z_23{ii},k_z{ii},w_z_23{ii}] = g2_cumulative_avg(out_corrs{ii}.g23,0.1);
    [k_z_34{ii},k_z{ii},w_z_34{ii}] = g2_cumulative_avg(out_corrs{ii}.g34,0.1);
    [k_z_14{ii},k_z{ii},w_z_14{ii}] = g2_cumulative_avg(out_corrs{ii}.g14,0.1);
end

%%
for ii = 1:length(phi_vec)      
    %% Separate data into the four ports
    ports = {};
    [ports.top_left, ports.top_right] = separate_ports(out_data{ii}.top_halo,0);
    [ports.bottom_left, ports.bottom_right] = separate_ports(out_data{ii}.bottom_halo,0);
    
    %% Quantum correlator E
    opts_E.calc_err = true;
    opts_E.plots = false;
    opts_E.verbose = false;
    opts_E.fit = false;
    opts_E.norm = false; %use normalised or unnormalised data
%     opts_E.sample_proportion = 1.0;%0.1;
    lambda = 0.1;
%     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
    opts_E.delta_kd = [4.2e-3.*4,1e-3,6e-3.*4];%
    opts_E.dim = 2;
    opts_E.sample_proportion = 1.0;
    opts_E.num_samp_rep = 50;
    
    %% BACK TO BACK (in the same halo)
    corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    
    global_sample_portion = 1.0;
    
    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    dkr = (dkx.*dky.*dkz).^(1/3);
    
    % variables for calculating the error
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.type='1d_cart_bb';%'radial_bb';%
    corr_opts.plots = true;
    corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
    corr_opts.fit = false;
    corr_opts.calc_err = do_g2_err;
    corr_opts.one_d_dimension = 2;
    corr_opts.two_d_dimensions = [2,3];
    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
    
    one_d_range=0.05;

    num_pts_cart = round(one_d_range./opts_E.delta_kd(corr_opts.one_d_dimension));
    num_pts_rad = round(one_d_range./dkr);
    num_pts_1 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(1)));
    num_pts_2 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(2)));
    corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
    corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);
    corr_opts.two_d_edges = {linspace(-one_d_range,one_d_range,num_pts_1)',linspace(-one_d_range,one_d_range,num_pts_2)'};
    corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';
    
    corr_opts.rad_smoothing=nan;
    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;
    
    corr_opts.norm_samp_factor=1500;%1500;
    corr_opts.sample_proportion=1.0;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    
    corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
    
    %% TOP HALO BACK TO BACK

        
        [E_val(ii), corrs.ports] = E(ports,opts_E);
        
        out_corrs{ii} = corrs.ports;
end

for ii = 1:num_phi
    [k_x_12{ii},k_x{ii},w_x_12{ii}] = g2_cumulative_avg(out_corrs{ii}.g12,0.1);
    [k_x_23{ii},k_x{ii},w_x_23{ii}] = g2_cumulative_avg(out_corrs{ii}.g23,0.1);
    [k_x_34{ii},k_x{ii},w_x_34{ii}] = g2_cumulative_avg(out_corrs{ii}.g34,0.1);
    [k_x_14{ii},k_x{ii},w_x_14{ii}] = g2_cumulative_avg(out_corrs{ii}.g14,0.1);
end

%%


for ii = 1:length(phi_vec)      
    %% Separate data into the four ports
    ports = {};
    [ports.top_left, ports.top_right] = separate_ports(out_data{ii}.top_halo,0);
    [ports.bottom_left, ports.bottom_right] = separate_ports(out_data{ii}.bottom_halo,0);
    
    %% Quantum correlator E
    opts_E.calc_err = true;
    opts_E.plots = false;
    opts_E.verbose = false;
    opts_E.fit = false;
    opts_E.norm = false; %use normalised or unnormalised data
%     opts_E.sample_proportion = 1.0;%0.1;
    lambda = 0.1;
%     opts_E.delta_kd = [4.2e-3,1.5e-3,6e-3].*4.*lambda;%[4e-3,1.3e-3,6e-3].*1;%[4e-3,1.3e-3,6e-3].*2;%[5e-3,2.*3e-3,10.*3e-3];%[3e-3,3e-3,3e-3];% volume widths in dimensions z x y used to calculate correlations
    opts_E.delta_kd = [4.2e-3.*4,1.5e-3.*4,1e-3];%
    opts_E.dim = 3;
    opts_E.sample_proportion = 1.0;
    opts_E.num_samp_rep = 50;
    
    %% BACK TO BACK (in the same halo)
    corr_opts.verbose = false;
    corr_opts.print_update = false;
    corr_opts.timer=false;
    
    global_sample_portion = 1.0;
    
    dkx = opts_E.delta_kd(2);
    dky = opts_E.delta_kd(3);
    dkz = opts_E.delta_kd(1);
    dkr = (dkx.*dky.*dkz).^(1/3);
    
    % variables for calculating the error
    corr_opts.samp_frac_lims=[0.65,0.9];
    corr_opts.num_samp_frac=5;
    corr_opts.num_samp_rep=5;
    
    corr_opts.attenuate_counts=1;
    corr_opts.type='1d_cart_bb';%'radial_bb';%
    corr_opts.plots = true;
    corr_opts.fig=['top halo bb corr ',num2str(phi_vec(ii))];
    corr_opts.fit = false;
    corr_opts.calc_err = do_g2_err;
    corr_opts.one_d_dimension = 2;
    corr_opts.two_d_dimensions = [2,3];
    corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];
    
    one_d_range=0.05;

    num_pts_cart = round(one_d_range./opts_E.delta_kd(corr_opts.one_d_dimension));
    num_pts_rad = round(one_d_range./dkr);
    num_pts_1 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(1)));
    num_pts_2 = round(one_d_range./opts_E.delta_kd(corr_opts.two_d_dimensions(2)));
    corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
    corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);
    corr_opts.two_d_edges = {linspace(-one_d_range,one_d_range,num_pts_1)',linspace(-one_d_range,one_d_range,num_pts_2)'};
    corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';
    
    corr_opts.rad_smoothing=nan;
    corr_opts.direction_labels = {'z','x','y'};
    corr_opts.low_mem=true;
    
    corr_opts.norm_samp_factor=1500;%1500;
    corr_opts.sample_proportion=1.0;%1.0;%0.65;%1500;
    corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')
    corr_opts.do_pre_mask=false;
    corr_opts.sorted_dir=nan;
    corr_opts.sort_norm=0;
    
    corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
    
    %% TOP HALO BACK TO BACK

        
        [E_val(ii), corrs.ports] = E(ports,opts_E);
        
        out_corrs{ii} = corrs.ports;
end

for ii = 1:num_phi
    [k_y_12{ii},k_y{ii},w_y_12{ii}] = g2_cumulative_avg(out_corrs{ii}.g12,0.1);
    [k_y_23{ii},k_y{ii},w_y_23{ii}] = g2_cumulative_avg(out_corrs{ii}.g23,0.1);
    [k_y_34{ii},k_y{ii},w_y_34{ii}] = g2_cumulative_avg(out_corrs{ii}.g34,0.1);
    [k_y_14{ii},k_y{ii},w_y_14{ii}] = g2_cumulative_avg(out_corrs{ii}.g14,0.1);
end
dt = '20210824';
save(['g2amp_vs_dk_',dt,'.mat'],'k_z_12','k_z_23','k_z_34','k_z_14','k_x_12','k_x_23','k_x_34','k_x_14','k_y_12','k_y_23','k_y_34','k_y_14','w_x_12','w_x_23','w_x_34','w_x_14','w_y_12','w_y_23','w_y_34','w_y_14','w_z_12','w_z_23','w_z_34','w_z_14','k_x','k_y','k_z')