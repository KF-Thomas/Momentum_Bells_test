function [top_halo_corr_bb,btm_halo_corr_bb,btw_halo_corr_bb,corr_amp] = g2_around_halo(top_halo,bottom_halo,global_opts)
%g2 around halo
%% general settings
corr_opts.verbose = false;
corr_opts.print_update = false;
corr_opts.timer=false;

corr_opts.plots = global_opts.plots;
corr_opts.fit = global_opts.fit;
corr_opts.calc_err = global_opts.calc_err;

global_sample_portion = 0.04;%1.0;%0.05;%0.08;%0.5;%1.0;%

% variables for calculating the error
corr_opts.samp_frac_lims=[0.65,0.9];
corr_opts.num_samp_frac=5;
corr_opts.num_samp_rep=5;

corr_opts.attenuate_counts=1;

%volume widths
dkx = global_opts.delta_kd(2);
dky = global_opts.delta_kd(3);
dkz = global_opts.delta_kd(1);
dkr = (dkx.*dky.*dkz).^(1/3);

%% BACK TO BACK (in the same halo)
%chose method of correlation calculation
corr_opts.type='1d_cart_bb';%'1d_vol_bb';%'2d_cart_bb';%'radial_bb';%
corr_density='one_d_corr_density';
corr_opts.bin_lims = 6;
corr_opts.one_d_dimension = 3;
corr_opts.two_d_dimensions = [2,3];
corr_opts.do_norm = true;

corr_opts.one_d_window=[[-1,1].*dkz;[-1,1].*dkx;[-1,1].*dky];

%setup grid spacing
one_d_range=0.05;%0.01;%0.01;%0.017;%0.09;%0.075;%0.02%0.03
% one_d_range=0.16;

num_pts_cart = round(one_d_range./global_opts.delta_kd(corr_opts.one_d_dimension));
num_pts_rad = round(one_d_range./dkr);
num_pts_1 = round(one_d_range./global_opts.delta_kd(corr_opts.two_d_dimensions(1)));
num_pts_2 = round(one_d_range./global_opts.delta_kd(corr_opts.two_d_dimensions(2)));
corr_opts.redges=sqrt(linspace(0^2,one_d_range^2,num_pts_rad));
corr_opts.one_d_edges = linspace(-one_d_range,one_d_range,num_pts_cart.*2);
corr_opts.two_d_edges = {linspace(-one_d_range,one_d_range,num_pts_1)',linspace(-one_d_range,one_d_range,num_pts_2)'};
corr_opts.edges=linspace(-1,1)';%corr_opts.edges=linspace(-1,-0.8)';

corr_opts.rad_smoothing=nan;%0.0003;%

corr_opts.direction_labels = {'z','x','y'};
corr_opts.low_mem=true;

%normalisation settings
corr_opts.norm_samp_factor=1500;%1500;
corr_opts.sample_proportion=global_sample_portion;%0.65;%1500;
corr_opts.sampling_method='complete';%'basic';%method for sampling uncorrelated pairs (either 'basic' or 'complete')

corr_opts.do_pre_mask=false;
corr_opts.sorted_dir=nan;
corr_opts.sort_norm=0;

corr_opts.gaussian_fit = true; %ensure it always uses a gaussian fit
corr_opts.param_num = 4;

%% bin in theta and phi
num_bins_theta = 24;
num_bins_phi = 9;
ang_lim = 15*pi/180;
theta_vec = linspace(0,pi,num_bins_theta);
phi_vec = linspace(pi/5,4/5*pi,num_bins_phi);%linspace(pi/6,5/6*pi,num_bins_phi);
% cart2sph(temp(:,1),temp(:,2),temp(:,3))
vec_norm = @(x) sqrt(sum(x.^2,2));
% counts = top_halo.counts_vel;
for kk = 1:num_bins_theta
    for jj = 1:num_bins_phi
        for ii = 1:length(top_halo.counts_txy)
            this_zxy_top = top_halo.counts_vel{ii};
            this_zxy_btm = bottom_halo.counts_vel{ii};
            % this_zxy = [1,1,4;1,1,1];
            theta = theta_vec(kk);
            phi = phi_vec(jj);
            norm_vec = [cos(phi),cos(theta).*sin(phi),sin(theta).*sin(phi)];
            u = repmat(norm_vec,[size(this_zxy_top,1),1]);
            current_angles_top = acos(dot(u',this_zxy_top')'./(vec_norm(u).*vec_norm(this_zxy_top)));
            ub = repmat(norm_vec,[size(this_zxy_btm,1),1]);
            current_angles_btm = acos(dot(ub',this_zxy_btm')'./(vec_norm(ub).*vec_norm(this_zxy_btm)));
            
            
            current_mask_top = (current_angles_top<ang_lim) | (current_angles_top>(pi-ang_lim));
            current_mask_btm = (current_angles_btm<ang_lim) | (current_angles_btm>(pi-ang_lim));
            top_counts{ii} = this_zxy_top(current_mask_top,:);
            btm_counts{ii} = this_zxy_btm(current_mask_btm,:);
        end
        %run correlator on each bin
        top_halo_corr_bb{kk,jj}=calc_any_g2_type(corr_opts,top_counts);
        btm_halo_corr_bb{kk,jj}=calc_any_g2_type(corr_opts,btm_counts);
        
        both_counts = [top_counts;btm_counts];
        btw_halo_corr_bb{kk,jj}=calc_any_g2_type(corr_opts,both_counts);
        
        mid_pt = ceil(length(top_halo_corr_bb{kk,jj}.norm_g2.g2_amp)/2);
        g2_amp_top(kk,jj) = top_halo_corr_bb{kk,jj}.norm_g2.g2_amp(mid_pt);
        g2_amp_btm(kk,jj) = btm_halo_corr_bb{kk,jj}.norm_g2.g2_amp(mid_pt);
        g2_amp_btw(kk,jj) = btw_halo_corr_bb{kk,jj}.norm_g2.g2_amp(mid_pt);
        
        G2_amp_top(kk,jj) = top_halo_corr_bb{kk,jj}.in_shot_corr.(corr_density)(mid_pt);
        G2_amp_btm(kk,jj) = btm_halo_corr_bb{kk,jj}.in_shot_corr.(corr_density)(mid_pt);
        G2_amp_btw(kk,jj) = btw_halo_corr_bb{kk,jj}.in_shot_corr.(corr_density)(mid_pt);
          
    end
end
 
        corr_amp.g2_amp_top = g2_amp_top;
        corr_amp.g2_amp_btm = g2_amp_btm;
        corr_amp.g2_amp_btw = g2_amp_btw;
        
        corr_amp.G2_amp_top = G2_amp_top;
        corr_amp.G2_amp_btm = G2_amp_btm;
        corr_amp.G2_amp_btw = G2_amp_btw;
        
        E_sph = (corr_amp.G2_amp_top+corr_amp.G2_amp_btm-2.*corr_amp.G2_amp_btw)./(corr_amp.G2_amp_top+corr_amp.G2_amp_btm+2.*corr_amp.G2_amp_btw);
        
        corr_amp.E_sph = E_sph;

if global_opts.global_plots
    stfig('g2 amp top halo');
    surf(theta_vec,phi_vec,g2_amp_top')
    xlabel('$\theta$')
    ylabel('$\phi$')
    zlabel('$g^{(2)}_{BB}(0)$')
    stfig('g2 amp bottom halo');
    surf(theta_vec,phi_vec,g2_amp_btm')
    xlabel('$\theta$')
    ylabel('$\phi$')
    zlabel('$g^{(2)}_{BB}(0)$')
    stfig('g2 amp between halos');
    surf(theta_vec,phi_vec,g2_amp_btw')
    xlabel('$\theta$')
    ylabel('$\phi$')
    zlabel('$g^{(2)}_{BB}(0)$')
    stfig('E amp around halo');
    surf(theta_vec,phi_vec,E_sph.')
    xlabel('$\theta$')
    ylabel('$\phi$')
    zlabel('$E(\Phi)$')
end
end



