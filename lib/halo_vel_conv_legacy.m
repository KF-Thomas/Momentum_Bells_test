function out_halo = halo_vel_conv(data,opts_vel_conv)
d = -opts_vel_conv.const.fall_distance;
g0 = -opts_vel_conv.const.g0;
tf = sqrt(2*-d/-g0);%arrival time of zero velocity particles
if opts_vel_conv.visual
    stfig(opts_vel_conv.title);
    clf
    xlabel('\(v_x\)')
    ylabel('\(v_y\)')
    zlabel('\(v_z\)')
    hold on
    axis equal
end
num_shots = length(data.shot_num);
for this_idx = 1:num_shots % Loop over all shots
    this_txy = data.counts_txy{this_idx};
    % center the data such that the BEC are at antipodal points
    this_centre = ones(1,3);
    dt = (-opts_vel_conv.bec_center.north(this_idx, 1)+opts_vel_conv.bec_center.south(this_idx, 1));
    new_t = dt/2+1/2.*sqrt(dt^2+8*d/g0);
    this_centre(1) = opts_vel_conv.bec_center.south(this_idx, 1)-new_t;
    dxy = -opts_vel_conv.bec_center.north(this_idx, 2:3)+opts_vel_conv.bec_center.south(this_idx, 2:3);
    new_xy = -(new_t-dt)/(2.*new_t-dt).*dxy;
    this_centre(2:3) = opts_vel_conv.bec_center.north(this_idx, 2:3)-new_xy;
    
    centred_counts = this_txy - this_centre;
    txy_top = opts_vel_conv.bec_center.north(this_idx, :)- this_centre;
    txy_mid = opts_vel_conv.bec_center.south(this_idx, :)- this_centre;
    txy_bec_top = [txy_top-opts_vel_conv.bec_width.north(this_idx, :);txy_top+opts_vel_conv.bec_width.north(this_idx, :)];
    txy_bec_mid = [txy_mid-opts_vel_conv.bec_width.south(this_idx, :);txy_mid+opts_vel_conv.bec_width.south(this_idx, :)];

    % Convert to kspace
    this_outtime = 0;%- 0.418707;%this_centre(1)
    v_top = txy_to_vel(txy_top, this_outtime, opts_vel_conv.const.g0, opts_vel_conv.const.fall_distance);
    v_mid = txy_to_vel(txy_mid, this_outtime, opts_vel_conv.const.g0, opts_vel_conv.const.fall_distance);
    v_bec_top = txy_to_vel(txy_bec_top, this_outtime, opts_vel_conv.const.g0, opts_vel_conv.const.fall_distance);
    v_bec_mid = txy_to_vel(txy_bec_mid, this_outtime, opts_vel_conv.const.g0, opts_vel_conv.const.fall_distance);
    v_radius = norm(v_top-v_mid)./2;
    v_rxy = sqrt(v_top(2).^2+v_top(3).^2);%radius in the vx vy dims
    theta = acos(v_top(3)./v_rxy)*180/pi;
    if isnan(theta)
        theta=0;
    end
    phi = asin(v_rxy/v_radius)*180/pi;
    phix = asin(v_top(2)/v_radius)*180/pi;
    phiy = asin(v_top(3)/v_radius)*180/pi;
    v_top_rot = v_top*rotx(theta)'*roty(-phi)';
    v_mid_rot = v_mid*rotx(theta)'*roty(-phi)';
    if v_radius>opts_vel_conv.v_thresh
        continue
    end
    v_zxy = txy_to_vel(centred_counts, this_outtime, opts_vel_conv.const.g0, opts_vel_conv.const.fall_distance);
%     v_zxy = v_zxy*rotx(theta)'*roty(-phi)';%rotate the BEC to the north and south poles
    v_zxy = v_zxy*roty(-phi)';%rotate the BEC to the north and south poles
    % mask radial
    radius_mask = (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)<(v_radius.*opts_vel_conv.v_mask(2)).^2 ...
    & (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)>(v_radius.*opts_vel_conv.v_mask(1)).^2;
    v_zxy = v_zxy(radius_mask,:);
    %mask out the BEC's
    v_masked = mask_square(v_zxy,v_bec_top',1);
    v_masked = mask_square(v_masked,v_bec_mid',1);
    %z mask
    z_lim = [opts_vel_conv.z_mask;-inf,inf;-inf,inf].*v_radius;
    v_masked = mask_square(v_masked,z_lim,0);
    %add the data to the structure
    out_halo.counts_txy{this_idx} = this_txy;
    out_halo.num_counts(this_idx) = size(v_masked,1);%size(this_txy,1);
    out_halo.counts_vel{this_idx} = v_masked;
    out_halo.counts_vel_unmasked{this_idx} = v_zxy;
    out_halo.rad(this_idx) = v_radius;
    out_halo.ang(this_idx,:) = [theta,phi];
    out_halo.bec_pos(this_idx,:) = [v_top,v_mid];
    out_halo.counts_vel_norm{this_idx} = v_masked./v_radius;
    %mask out the BEC's
    if opts_vel_conv.visual && rand(1)<opts_vel_conv.plot_percentage
        scatter3(v_masked(:,2),v_masked(:,3),v_masked(:,1),'k.')
    end
end
out_halo.shot_num = data.shot_num;
end