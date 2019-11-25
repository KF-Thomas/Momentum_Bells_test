function plot_checks(halos,plot_opts)
    top_halo = halos.top_halo;
    bottom_halo = halos.bottom_halo;
    bec_masked = halos.bec;
    
    v_top_zxy = cell2mat(top_halo.counts_vel_norm);
    r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
    N_top = top_halo.num_counts;
    
    v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm);
    r_dist_btm = sqrt(v_btm_zxy(:,1).^2+v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2);
    N_btm = bottom_halo.num_counts;
    
    stfig('halo comparison')
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
    plot(r_hist_top.bin.centers,r_hist_top.counts.smooth,'linewidth',1.5)
    hold on
    plot(r_hist_btm.bin.centers,r_hist_btm.counts.smooth,'linewidth',1.5)
    xlabel('r')
    ylabel('Freq')
    xlim([min([r_hist_top.bin.centers;r_hist_btm.bin.centers]),...
        max([r_hist_top.bin.centers;r_hist_btm.bin.centers])])
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
    stfig('density of top halo');
    subplot(1,3,1)
    ndhist(v_top_zxy(:,2:3));
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    subplot(1,3,2)
    ndhist(v_top_zxy(:,[1,3]));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_y$')
    subplot(1,3,3)
    ndhist(v_top_zxy(:,1:2));
    axis equal
    xlabel('$v_z$')
    ylabel('$v_x$')
end