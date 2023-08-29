function out_halo = halo_vel_conv(data,opts_vel_conv)
d = opts_vel_conv.const.fall_distance;
g0 = opts_vel_conv.const.g0;
tf = sqrt(2*d/g0);%arrival time of zero velocity particles
if opts_vel_conv.visual
    stfig(opts_vel_conv.title);
    clf
    xlabel('\(v_x\)')
    ylabel('\(v_y\)')
    zlabel('\(v_z\)')
    hold on
    axis equal
end

if ~isfield(opts_vel_conv,'phi_correction')
    opts_vel_conv.phi_correction = [0 0];
end

if isfield(opts_vel_conv,'ang_lim')
    ang_lim = opts_vel_conv.ang_lim;
else
    ang_lim = 60;
end

num_shots = length(data.shot_num);
for this_idx = 1:num_shots % Loop over all shots
    
    %     dt = (-opts_vel_conv.bec_center.north(this_idx, 1)+opts_vel_conv.bec_center.south(this_idx, 1));
    %     new_t = dt/2+1/2.*sqrt(dt^2+8*d/g0);
    %     this_centre(1) = opts_vel_conv.bec_center.south(this_idx, 1)-new_t;
    %     dxy = -opts_vel_conv.bec_center.north(this_idx, 2:3)+opts_vel_conv.bec_center.south(this_idx, 2:3);
    %     new_xy = -(new_t-dt)/(2.*new_t-dt).*dxy;
    %     this_centre(2:3) = opts_vel_conv.bec_center.north(this_idx, 2:3)-new_xy;
    %% INTIAL DATA
    %current txy data and zero momentum point
    this_centre = opts_vel_conv.center(this_idx,:); %center the zero momentum point (for the most part the top BEC)
    this_txy = data.counts_txy{this_idx};
    
    % read in the txy and centering data
    centred_counts = this_txy - this_centre;
    txy_north = opts_vel_conv.bec_center.north(this_idx, :)- this_centre;
    txy_south = opts_vel_conv.bec_center.south(this_idx, :)- this_centre;
    txy_bec_top = [txy_north-opts_vel_conv.bec_width.north(this_idx, :);txy_north+opts_vel_conv.bec_width.north(this_idx, :)];
    txy_bec_mid = [txy_south-opts_vel_conv.bec_width.south(this_idx, :);txy_south+opts_vel_conv.bec_width.south(this_idx, :)];
    
    %% VELOCITY CONVERSION
    % Convert the poles (BEC's) to kspace
    this_outtime = -tf;
    vel_shift = txy_to_vel(txy_north, this_outtime, g0, d);%+opts_vel_conv.centering_correction;
    v_north = txy_to_vel(txy_north, this_outtime, g0, d)-vel_shift;
    v_south = txy_to_vel(txy_south, this_outtime, g0, d)-vel_shift;
    v_bec_top = txy_to_vel(txy_bec_top, this_outtime, g0, d);
    v_bec_mid = txy_to_vel(txy_bec_mid, this_outtime, g0, d);
    v_radius = norm(v_north-v_south)./2;
    
    % Check if radius is a reasonable value
    if v_radius>opts_vel_conv.v_thresh
        continue
    end
    
    sum_v_z = v_north(1)+v_south(1); %sum of vz values
    z_sign = sign(sum_v_z);
    v_zx = sqrt(sum_v_z.^2+(v_north(2)+v_south(2)).^2);
    phix = atan((v_north(2)+v_south(2))/sum_v_z)*180/pi;
    phiy = atan((v_north(3)+v_south(3))/v_zx)*180/pi;
    
    % Positions in velocity space of the north and south poles of the halo
    v_north_rot = v_north*rotz(-phix)'*roty(-phiy)'- v_radius.*[z_sign,0,0];
    v_south_rot = v_south*rotz(-phix)'*roty(-phiy)'- v_radius.*[z_sign,0,0];
    
    time_mask = centred_counts(:,1)>this_outtime;
    centred_counts_masked = centred_counts(time_mask,:);
    % Convert raw counts to velocity space
    v_zxy = txy_to_vel(centred_counts_masked, this_outtime, g0, d)-vel_shift;
    v_zxy = v_zxy*rotz(-phix)'*roty(-phiy)';%rotate the BEC to the north and south poles
    v_zxy = v_zxy - [z_sign*v_radius,0,0];%*rotz(phix)'*roty(phiy)'; %shift into center of mass frame
    phix_correction = opts_vel_conv.phi_correction(1);
    phiy_correction = opts_vel_conv.phi_correction(2);
    v_zxy = (v_zxy - opts_vel_conv.centering_correction)*rotz(-phix_correction)'*roty(-phiy_correction)';
    
    v_zxy_unmasked = v_zxy;
    
    %% MASKING
    % mask radial
    radius_mask = (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)<(opts_vel_conv.v_mask(2)).^2 ...
        & (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)>(opts_vel_conv.v_mask(1)).^2;
    v_zxy = v_zxy(radius_mask,:);
    
    if isfield(opts_vel_conv,'vxy_mask')
    xy_rad_mask = (v_zxy(:,2).^2+v_zxy(:,3).^2)<(opts_vel_conv.vxy_mask(2)).^2 ...
        & (v_zxy(:,2).^2+v_zxy(:,3).^2)>(opts_vel_conv.vxy_mask(1)).^2;
    v_zxy = v_zxy(logical(xy_rad_mask),:);
    end

    if isfield(opts_vel_conv,'theta_mask')
        [theta,rho,z] = cart2pol(v_zxy(:,2),v_zxy(:,3),v_zxy(:,1));
        theta_mask = false(size(theta));
        for jj = 1:size(opts_vel_conv.theta_mask,1)
        theta_mask = theta_mask | (theta>opts_vel_conv.theta_mask(jj,1) & theta<opts_vel_conv.theta_mask(jj,2));
        end
        if isfield(opts_vel_conv,'theta_mask_polarity') && opts_vel_conv.theta_mask_polarity
        v_zxy = v_zxy(theta_mask,:);
        else
            v_zxy = v_zxy(~theta_mask,:);
        end
    end
    
    ang_mask = abs(180/pi*atan(v_zxy(:,1)./sqrt(v_zxy(:,2).^2+v_zxy(:,3).^2)))<ang_lim;
    v_zxy = v_zxy(ang_mask,:);
    
    %mask out the BEC's
    %     v_masked = mask_square(v_zxy,v_bec_top',1);
    %     v_masked = mask_square(v_masked,v_bec_mid',1);
    
    %z and y mask
    z_lim = [opts_vel_conv.z_mask;-inf,inf;-inf,inf].*v_radius;
    y_lim = [-inf,inf;-inf,inf;opts_vel_conv.y_mask].*v_radius;
    v_masked = mask_square(v_zxy,z_lim,0);
    v_masked = mask_square(v_masked,y_lim,0);
    
    %do some angular masking
    
    %% some data for post checking
    num_check = sum(this_txy(:,1)>1.5 & this_txy(:,1)<1.689898);
    
    top_halo_num = sum(this_txy(:,1)>1.76 & this_txy(:,1)<1.769);
    btm_halo_num = sum(this_txy(:,1)>1.747 & this_txy(:,1)<1.753);
    
    out_halo.top_halo_num(this_idx) = top_halo_num;
    out_halo.btm_halo_num(this_idx) = btm_halo_num;
    
    
    %% add the data to the structure
    out_halo.counts_txy{this_idx} = this_txy;
    out_halo.num_counts(this_idx) = size(v_masked,1);%size(this_txy,1);
    out_halo.counts_vel{this_idx} = v_masked;
    out_halo.counts_vel_unmasked{this_idx} = v_zxy;
    out_halo.rad(this_idx) = v_radius;
    out_halo.ang(this_idx,:) = [phix,phiy];
    out_halo.bec_pos(this_idx,:) = [v_north,v_south];
    out_halo.bec_pos_adj(this_idx,:) = [v_north_rot,v_south_rot];
    out_halo.counts_vel_norm{this_idx} = v_masked./v_radius;
    out_halo.v_zxy_unmasked{this_idx} = v_zxy_unmasked;
    out_halo.num_check(this_idx) = num_check;
    
    if opts_vel_conv.visual && rand(1)<opts_vel_conv.plot_percentage
        scatter3(v_masked(:,2),v_masked(:,3),v_masked(:,1),'k.')
    end
end
out_halo.shot_num = data.shot_num;
end