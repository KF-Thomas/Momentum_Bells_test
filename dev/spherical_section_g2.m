function out_corrs = spherical_section_g2(corr_opts,counts_vel)

%set up grid
theta_bins = corr_opts.theta_bins;
phi_bins = corr_opts.phi_bins;

kk = 1;
%split halos into bins
for ii = 1:size(counts_vel,2)
    for nn = 1:size(counts_vel,1)
        v_zxy = counts_vel{nn,ii};
        [theta,~] = cart2pol(v_zxy(:,2),v_zxy(:,3));
        phi = atan(v_zxy(:,1)./sqrt(v_zxy(:,2).^2+v_zxy(:,3).^2));

        phi_mask = (phi_bins(kk)<phi & phi<=phi_bins(kk+1)) |...
            ((-phi_bins(kk+1))<phi & phi<=(-phi_bins(kk)));

        theta_mask = false(size(theta));
        for jj = 1:size(theta_bins,1)
            theta_mask = theta_mask | (theta_bins(jj,1)<theta & theta<=theta_bins(jj,2))| ...
                ((theta_bins(jj,1)-pi)<theta & theta<=(theta_bins(jj,2)-pi));
        end

        if isfield(corr_opts,'theta_mask_polarity') && corr_opts.theta_mask_polarity == 0
            ang_mask = phi_mask & ~theta_mask;
        else
            ang_mask = phi_mask & theta_mask;
        end
        masked_counts_vel{nn,ii}  = v_zxy(ang_mask,:);
    end 
end

%calculate g2
out_corrs.data = calc_any_g2_type(corr_opts,masked_counts_vel);

end