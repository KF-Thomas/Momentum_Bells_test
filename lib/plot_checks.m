function plot_checks(halos,plot_opts)
d = plot_opts.const.fall_distance;
g0 = plot_opts.const.g0;
tf = sqrt(2*d/g0);%arrival time of zero velocity particles

    top_halo = halos.top_halo;
    bottom_halo = halos.bottom_halo;
    bec_masked = halos.bec;
    
    v_top_zxy = cell2mat(top_halo.counts_vel_norm);
    v_top_zxy_unnorm = cell2mat(top_halo.counts_vel);
    r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
    r_dist_top_unnorm = sqrt(v_top_zxy_unnorm(:,1).^2+v_top_zxy_unnorm(:,2).^2+v_top_zxy_unnorm(:,3).^2);
    N_top = top_halo.num_counts;
    
    v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm);
    v_btm_zxy_unnorm = cell2mat(bottom_halo.counts_vel);
    r_dist_btm = sqrt(v_btm_zxy(:,1).^2+v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2);
    r_dist_btm_unnorm = sqrt(v_btm_zxy_unnorm(:,1).^2+v_btm_zxy_unnorm(:,2).^2+v_btm_zxy_unnorm(:,3).^2);
    N_btm = bottom_halo.num_counts;
    
    stfig('Alignment of the BECs');
    clf
    subplot(5,1,1)
    plot(bec_masked.shot_num,bec_masked.centre_top(:,2:3))
    hold on
    plot(bec_masked.shot_num,bec_masked.centre_mid(:,2:3))
    plot(bec_masked.shot_num,bec_masked.centre_btm(:,2:3))
    ylabel('position (m)')
    xlabel('shot num')
    subplot(5,1,2)
    plot(bec_masked.shot_num,bec_masked.centre_top(:,1))
    hold on
    plot(bec_masked.shot_num,bec_masked.centre_mid(:,1))
    plot(bec_masked.shot_num,bec_masked.centre_btm(:,1))
    ylabel('position (t)')
    xlabel('shot num')
    subplot(5,1,3)
    scatter(bec_masked.centre_top(:,1),bec_masked.centre_top(:,3),'.')
    hold on
    scatter(bec_masked.centre_mid(:,1),bec_masked.centre_mid(:,3),'.')
    scatter(bec_masked.centre_btm(:,1),bec_masked.centre_btm(:,3),'.')
    scatter(bec_masked.centre_top(:,1),bec_masked.centre_top(:,2),'.')
    scatter(bec_masked.centre_mid(:,1),bec_masked.centre_mid(:,2),'.')
    scatter(bec_masked.centre_btm(:,1),bec_masked.centre_btm(:,2),'.')
    xlabel('arrival time (s)')
    ylabel('x/y position (m)')
    subplot(5,1,4)
    this_outtime = bec_masked.centre_mid(:,1) - tf;
    top_txy = [bec_masked.centre_top(:,1),bec_masked.centre_top(:,2)-bec_masked.centre_mid(:,2),...
        bec_masked.centre_top(:,3)-bec_masked.centre_mid(:,3)];
    mid_txy = [bec_masked.centre_mid(:,1),bec_masked.centre_mid(:,2)-bec_masked.centre_mid(:,2),...
        bec_masked.centre_mid(:,3)-bec_masked.centre_mid(:,3)];
    btm_txy = [bec_masked.centre_btm(:,1),bec_masked.centre_btm(:,2)-bec_masked.centre_mid(:,2),...
        bec_masked.centre_btm(:,3)-bec_masked.centre_mid(:,3)];
    top_v = txy_to_vel(top_txy, this_outtime, g0, d);
    mid_v = txy_to_vel(mid_txy, this_outtime, g0, d);
    btm_v = txy_to_vel(btm_txy, this_outtime, g0, d);
    scatter(mid_v(:,1),mid_v(:,3))
    hold on
    theta_top_y = atan(top_v(:,3)./top_v(:,1));
    theta_btm_y = atan(btm_v(:,3)./btm_v(:,1));
    theta_top_x = atan(top_v(:,2)./top_v(:,1));
    theta_btm_x = atan(btm_v(:,2)./btm_v(:,1));
%     rot_mat = 
    scatter(top_v(:,1),top_v(:,3))
    scatter(btm_v(:,1),btm_v(:,3))
    scatter(top_v(:,1),top_v(:,2))
    scatter(mid_v(:,1),mid_v(:,2))
    scatter(btm_v(:,1),btm_v(:,2))
%     scatter(btm_v(:,1),sqrt(btm_v(:,2).^2+btm_v(:,3).^2),'.')
%     scatter(top_v(:,1),sqrt(top_v(:,2).^2+top_v(:,3).^2),'.')
    xlabel('\(v_z\)')
    ylabel('\(v_{x/y}\)')
    subplot(5,2,9)
    histogram(theta_top_y.*180/pi,100)
    hold on
    histogram(theta_btm_y.*180/pi,100)
    fprintf('phi_{t,y} ~ %s\n',string_value_with_unc(mean(theta_top_y.*180/pi),std(theta_top_y.*180/pi),'type','b','separator',0))
    fprintf('phi_{b,y} ~ %s\n',string_value_with_unc(mean(theta_btm_y.*180/pi),std(theta_btm_y.*180/pi),'type','b','separator',0))
    xlabel('Rot ang (degree)')
    ylabel('Count')
    legend('\(\phi_{t,y}\)','\(\phi_{b,y}\)')
    subplot(5,2,10)
    histogram(theta_top_x.*180/pi,100)
    hold on
    histogram(theta_btm_x.*180/pi,100)
    fprintf('phi_{t,x} ~ %s\n',string_value_with_unc(mean(theta_top_x.*180/pi),std(theta_top_x.*180/pi),'type','b','separator',0))
    fprintf('phi_{b,x} ~ %s\n',string_value_with_unc(mean(theta_btm_x.*180/pi),std(theta_btm_x.*180/pi),'type','b','separator',0))
    xlabel('Rot ang (degree)')
    ylabel('Count')
    legend('\(\phi_{t,x}\)','\(\phi_{b,x}\)')
    
    stfig('distribution of halo radius');
    clf
    vr_hist_top=smooth_hist(top_halo.rad,'sigma',0.00003);
    vr_hist_btm=smooth_hist(bottom_halo.rad,'sigma',0.00003);
    plot(vr_hist_top.bin.centers,vr_hist_top.counts.smooth,'linewidth',1.5)
    hold on
    plot(vr_hist_btm.bin.centers,vr_hist_btm.counts.smooth,'linewidth',1.5)
    xlabel('radius of halo')
    ylabel('Freq')
    legend('top','btm')
    
    stfig('halo comparison');
    clf
    plot_mask_top = rand(size(v_top_zxy,1),1)<0.05;
    plot_mask_btm = rand(size(v_btm_zxy,1),1)<0.05;
    scatter3(v_top_zxy(plot_mask_top,2),v_top_zxy(plot_mask_top,3),v_top_zxy(plot_mask_top,1),'r.')
    hold on
    scatter3(v_btm_zxy(plot_mask_btm,2),v_btm_zxy(plot_mask_btm,3),v_btm_zxy(plot_mask_btm,1),'b.')
    axis equal
    
    
    stfig('radial distribution');
    clf
    r_hist_top=smooth_hist(r_dist_top,'sigma',0.0001);
    r_hist_btm=smooth_hist(r_dist_btm,'sigma',0.0001);
    r_hist_top_un=smooth_hist(r_dist_top_unnorm,'sigma',0.0001);
    r_hist_btm_un=smooth_hist(r_dist_btm_unnorm,'sigma',0.0001);
    subplot(2,1,1)
    plot(r_hist_top.bin.centers,r_hist_top.counts.smooth,'linewidth',1.5)
    hold on
    plot(r_hist_btm.bin.centers,r_hist_btm.counts.smooth,'linewidth',1.5)
    xlabel('r')
    ylabel('Freq')
    xlim([min([r_hist_top.bin.centers;r_hist_btm.bin.centers]),...
        max([r_hist_top.bin.centers;r_hist_btm.bin.centers])])
    legend('top','btm')
    subplot(2,1,2)
    plot(r_hist_top_un.bin.centers,r_hist_top_un.counts.smooth,'linewidth',1.5)
    hold on
    plot(r_hist_btm_un.bin.centers,r_hist_btm_un.counts.smooth,'linewidth',1.5)
    xlabel('r')
    ylabel('Freq')
    xlim([min([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers]),...
        max([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers])])
    legend('top','btm')
    stfig('Counts in halo distribution');
    clf
    hold on
    N_hist_top=smooth_hist(N_top,'lims',[0,100],'sigma',1);
    N_hist_btm=smooth_hist(N_btm,'lims',[0,100],'sigma',1);
    Nplot1=plot(N_hist_top.bin.centers,N_hist_top.counts.smooth,'linewidth',1.5);
    Nplot2=plot(N_hist_btm.bin.centers,N_hist_btm.counts.smooth,'linewidth',1.5);
    xlim([min([N_hist_top.bin.centers;N_hist_btm.bin.centers]),...
        max([N_hist_top.bin.centers;N_hist_btm.bin.centers])])
    box on
    xlabel('N')
    ylabel('Freq')
    legend('top','btm')
    
    stfig('Halo Num History');
    plot(top_halo.shot_num,...
        top_halo.num_counts,...
        'kx-','LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('N in Halo')
    
    trans_frac = [bec_masked.trans_top';bec_masked.trans_mid';bec_masked.trans_btm';bec_masked.trans_oth'];
    stfig('Transfer Fraction History');
    plot(top_halo.shot_num,...
        trans_frac',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('Transfer Fraction')
    legend('top','mid','btm','oth')
    
    stfig('Vel Width History');
    subplot(3,1,1)
    plot(top_halo.shot_num,...
        bec_masked.vel_width_top','.',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('Velocity width (m/s)')
    title('Top')
    legend('wz','wx','wy')
    subplot(3,1,2)
    plot(top_halo.shot_num,...
        bec_masked.vel_width_mid','.',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('Velocity width (m/s)')
    title('Midle')
    legend('wz','wx','wy')
    subplot(3,1,3)
    plot(top_halo.shot_num,...
        bec_masked.vel_width_btm','.',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('Velocity width (m/s)')
    title('Bottom')
    legend('wz','wx','wy')
    
    stfig('Mode Occupancy and Correlation Amplitude');
    subplot(2,1,1)
    plot(top_halo.shot_num,...
        top_halo.m,'kx-',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('Mode Occupancy')
    subplot(2,1,2)
    plot(top_halo.shot_num,...
        top_halo.g2,'kx-',...
        'LineWidth',1.5)
    grid on
    h=gca;
    grid on    % turn on major grid lines
    grid minor % turn on minor grid lines
    % Set limits and grid spacing separately for the two directions:
    % Must set major grid line properties for both directions simultaneously:
    h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
    h.GridAlpha=1;  % the default is partially transparent
    h.GridColor=[0,0,0]; % here's the color for the major grid lines
    % Idem for minor grid line properties:
    h.MinorGridLineStyle='-';
    h.MinorGridAlpha=0.1;
    h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
    xlabel('Shot Number')
    ylabel('$g^2$')
 
    %getter better look at the density distribution
%     v_top_zxy = cell2mat(top_halo.counts_vel_unmasked);
    t=linspace(0,2*pi,1e4);
    stfig('density of top halo');
    clf
    subplot(2,3,1)
    ndhist(v_top_zxy(:,2:3));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,2)
    ndhist(v_top_zxy(:,[1,3]));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,3)
    ndhist(v_top_zxy(:,1:2));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
    
    subplot(2,3,4)
    ndhist(v_top_zxy_unnorm(:,2:3));
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(2,3,5)
    ndhist(v_top_zxy_unnorm(:,[1,3]));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(2,3,6)
    ndhist(v_top_zxy_unnorm(:,1:2));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
    
    stfig('density of bottom halo');
    clf
    subplot(1,3,1)
    ndhist(v_btm_zxy(:,2:3));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(1,3,2)
    ndhist(v_btm_zxy(:,[1,3]));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(1,3,3)
    ndhist(v_btm_zxy(:,1:2));
    hold on
    plot(cos(t),sin(t),'k')
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
end