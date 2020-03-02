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
num_shots = length(data.shot_num);
for this_idx = 1:num_shots % Loop over all shots
    this_txy = data.counts_txy{this_idx};
    % center the data such that the BEC are at antipodal points
    this_centre = opts_vel_conv.center(this_idx,:); %center the zero momentum point (for the most part the middle BEC)
%     dt = (-opts_vel_conv.bec_center.north(this_idx, 1)+opts_vel_conv.bec_center.south(this_idx, 1));
%     new_t = dt/2+1/2.*sqrt(dt^2+8*d/g0);
%     this_centre(1) = opts_vel_conv.bec_center.south(this_idx, 1)-new_t;
%     dxy = -opts_vel_conv.bec_center.north(this_idx, 2:3)+opts_vel_conv.bec_center.south(this_idx, 2:3);
%     new_xy = -(new_t-dt)/(2.*new_t-dt).*dxy;
%     this_centre(2:3) = opts_vel_conv.bec_center.north(this_idx, 2:3)-new_xy;
    
    centred_counts = this_txy - this_centre;
    txy_north = opts_vel_conv.bec_center.north(this_idx, :)- this_centre;
    txy_south = opts_vel_conv.bec_center.south(this_idx, :)- this_centre;
    txy_bec_top = [txy_north-opts_vel_conv.bec_width.north(this_idx, :);txy_north+opts_vel_conv.bec_width.north(this_idx, :)];
    txy_bec_mid = [txy_south-opts_vel_conv.bec_width.south(this_idx, :);txy_south+opts_vel_conv.bec_width.south(this_idx, :)];

    % Convert to kspace
    this_outtime = -tf;
    v_north = txy_to_vel(txy_north, this_outtime, g0, d);
    v_south = txy_to_vel(txy_south, this_outtime, g0, d);
    v_bec_top = txy_to_vel(txy_bec_top, this_outtime, g0, d);
    v_bec_mid = txy_to_vel(txy_bec_mid, this_outtime, g0, d);
    v_radius = norm(v_north-v_south)./2;
    v_rxy = sqrt(v_north(2).^2+v_north(3).^2);%radius in the vx vy dims
%     theta = acos(v_north(3)./v_rxy)*180/pi;
%     phi = asin(v_rxy/v_radius)*180/pi;
    v_z = v_north(1)+v_south(1);
    z_sign = sign(v_z);
    v_zx = sqrt(v_z.^2+(v_north(2)+v_south(2)).^2);
    phix = atan((v_north(2)+v_south(2))/v_z)*180/pi;
    phiy = atan((v_north(3)+v_south(3))/v_zx)*180/pi;
    try
    v_north_rot = v_north*rotz(-phix)'*roty(-phiy)'- v_radius.*[z_sign,0,0];
    v_south_rot = v_south*rotz(-phix)'*roty(-phiy)'- v_radius.*[z_sign,0,0];
    catch
        dum=0;
    end
    if v_radius>opts_vel_conv.v_thresh
        continue
    end
    v_zxy = txy_to_vel(centred_counts, this_outtime, g0, d);
    v_zxy = v_zxy*rotz(-phix)'*roty(-phiy)';%rotate the BEC to the north and south poles
    v_zxy(:,1) = v_zxy(:,1) - z_sign*v_radius;
    % mask radial
    radius_mask = (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)<(v_radius.*opts_vel_conv.v_mask(2)).^2 ...
    & (v_zxy(:,1).^2+v_zxy(:,2).^2+v_zxy(:,3).^2)>(v_radius.*opts_vel_conv.v_mask(1)).^2;
    v_zxy = v_zxy(radius_mask,:);
    %mask out the BEC's
%     v_masked = mask_square(v_zxy,v_bec_top',1);
%     v_masked = mask_square(v_masked,v_bec_mid',1);
    %z mask
    z_lim = [opts_vel_conv.z_mask;-inf,inf;-inf,inf].*v_radius;
    v_masked = mask_square(v_zxy,z_lim,0);
    %add the data to the structure
    out_halo.counts_txy{this_idx} = this_txy;
    out_halo.num_counts(this_idx) = size(v_masked,1);%size(this_txy,1);
    out_halo.counts_vel{this_idx} = v_masked;
    out_halo.counts_vel_unmasked{this_idx} = v_zxy;
    out_halo.rad(this_idx) = v_radius;
    out_halo.ang(this_idx,:) = [phix,phiy];
    out_halo.bec_pos(this_idx,:) = [v_north,v_south];
    out_halo.bec_pos_adj(this_idx,:) = [v_north_rot,v_south_rot];
    out_halo.counts_vel_norm{this_idx} = v_masked./v_radius;
    %mask out the BEC's
    if opts_vel_conv.visual && rand(1)<opts_vel_conv.plot_percentage
        scatter3(v_masked(:,2),v_masked(:,3),v_masked(:,1),'k.')
    end
end
out_halo.shot_num = data.shot_num;
end