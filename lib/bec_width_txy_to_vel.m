function bec = bec_width_txy_to_vel(bec,opts)
%converts the BEC from txy to velocity
bec.vel_width_top(:,1) = bec.width_top(:,1).*opts.g0;
bec.vel_width_top(:,2) = bec.width_top(:,2)./opts.fall_time;
bec.vel_width_top(:,3) = bec.width_top(:,3)./opts.fall_time;

bec.vel_width_mid(:,1) = bec.width_mid(:,1).*opts.g0;
bec.vel_width_mid(:,2) = bec.width_mid(:,2)./opts.fall_time;
bec.vel_width_mid(:,3) = bec.width_mid(:,3)./opts.fall_time;

bec.vel_width_btm(:,1) = bec.width_btm(:,1).*opts.g0;
bec.vel_width_btm(:,2) = bec.width_btm(:,2)./opts.fall_time;
bec.vel_width_btm(:,3) = bec.width_btm(:,3)./opts.fall_time;
end