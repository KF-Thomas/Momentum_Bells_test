function out_corrs = spherical_g2(corr_opts,counts_vel)

%set up grid
theta_bins = linspace(0,pi,corr_opts.num_theta_bins+1);
phi_bin_size = corr_opts.phi_max/(corr_opts.num_phi_bins-1/2);
phi_bins = [-(phi_bin_size/2),(phi_bin_size/2):phi_bin_size:corr_opts.phi_max];

nbins_ph = length(phi_bins);

%split halos into bins
for ii = 1:size(counts_vel,2)
    for nn = 1:size(counts_vel,1)
    v_zxy = counts_vel{nn,ii};
    [theta,rxy] = cart2pol(v_zxy(:,2),v_zxy(:,3));
    phi = atan(v_zxy(:,1)./sqrt(v_zxy(:,2).^2+v_zxy(:,3).^2));

    for jj = 1:(corr_opts.num_theta_bins)
        for kk = 1:(nbins_ph-1)
            if kk >1
            ang_mask = (theta_bins(jj)<theta & theta<=theta_bins(jj+1) ...
                & phi_bins(kk)<phi & phi<=phi_bins(kk+1)) | ...
                ((theta_bins(jj)-pi)<theta & theta<=(theta_bins(jj+1)-pi) ...
                & (-phi_bins(kk+1))<phi & phi<=(-phi_bins(kk)));
            else
             ang_mask = (theta_bins(jj)<theta & theta<=theta_bins(jj+1) ...
                & phi_bins(kk)<phi & phi<=phi_bins(kk+1)) | ...
                ((theta_bins(jj)-pi)<theta & theta<=(theta_bins(jj+1)-pi) ...
                & phi_bins(kk)<phi & phi<=phi_bins(kk+1));
            end

            %             area = abs(theta_bins(jj+1)-theta_bins(jj))*abs(sin(phi_bins(kk+1))-sin(phi_bins(kk)));

            %             sph_density(jj,kk) = sum(ang_mask)./area;
            masked_counts_vel{jj,kk}{nn,ii}  = v_zxy(ang_mask,:);

            phi_mean(kk) = mean(phi_bins(kk:(kk+1)));
        end
        theta_mean(jj) = mean(theta_bins(jj:(jj+1)));

    end
    end
end
out_corrs.theta_bins = theta_bins;
out_corrs.phi_bins = phi_bins;
out_corrs.theta_mid = theta_mean;
out_corrs.phi_mid = phi_mean;
%calculate g2
for jj = 1:(corr_opts.num_theta_bins)
    for kk = 1:(nbins_ph-1)
        out_corrs.data{jj,kk} = calc_any_g2_type(corr_opts,masked_counts_vel{jj,kk});
    end
end
end