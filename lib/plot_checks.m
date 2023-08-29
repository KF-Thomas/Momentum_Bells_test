function plot_checks(halos,plot_opts)
d = plot_opts.const.fall_distance;
g0 = plot_opts.const.g0;
tf = sqrt(2*d/g0);%arrival time of zero velocity particles

top_halo = halos.top_halo;
bottom_halo = halos.bottom_halo;
bec_masked = halos.bec;


%% find the widths of the halos in velocity space
top_halo.bec_vel_width = (mean(bec_masked.vel_width_top,2)+mean(bec_masked.vel_width_mid,2))./2;% add the average bec width
bottom_halo.bec_vel_width = (mean(bec_masked.vel_width_btm,2)+mean(bec_masked.vel_width_mid,2))./2;
%% find the mode number
opts.mode_num.qe = 0.08;
top_halo.m = halo_mode_occupancy(top_halo,opts.mode_num);
bottom_halo.m = halo_mode_occupancy(bottom_halo,opts.mode_num);
%% calculated expected correlation amplitude
top_halo.g2 = 2 + 1./top_halo.m;
bottom_halo.g2 = 2 + 1./bottom_halo.m;

v_top_zxy = cell2mat(top_halo.counts_vel_norm.');
v_top_zxy_unnorm = cell2mat(top_halo.counts_vel.');
r_dist_top = sqrt(v_top_zxy(:,1).^2+v_top_zxy(:,2).^2+v_top_zxy(:,3).^2);
r_dist_top_unnorm = sqrt(v_top_zxy_unnorm(:,1).^2+v_top_zxy_unnorm(:,2).^2+v_top_zxy_unnorm(:,3).^2);
N_top = top_halo.num_counts;

v_btm_zxy = cell2mat(bottom_halo.counts_vel_norm.');
v_btm_zxy_unnorm = cell2mat(bottom_halo.counts_vel.');
r_dist_btm = sqrt(v_btm_zxy(:,1).^2+v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2);
r_dist_btm_unnorm = sqrt(v_btm_zxy_unnorm(:,1).^2+v_btm_zxy_unnorm(:,2).^2+v_btm_zxy_unnorm(:,3).^2);
N_btm = bottom_halo.num_counts;

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
ylimit = max([r_hist_top_un.counts.smooth;r_hist_btm_un.counts.smooth]);
plot([0.130159/2 0.130159/2],[-0.1,ylimit.*2],'k-','linewidth',1.5)
ylim([0 ylimit.*1.1])
xlabel('r')
ylabel('Freq')
xlim([min([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers]),...
    max([r_hist_top_un.bin.centers;r_hist_btm_un.bin.centers])])
legend('top','btm','expected radius')
%     if plot_opts.only_dists
%         error('breaking out')
%     end
%
stfig('halo comparison');
clf
plot_mask_top = rand(size(v_top_zxy,1),1)<0.65;
plot_mask_btm = rand(size(v_btm_zxy,1),1)<0.65;
scatter3(v_top_zxy(plot_mask_top,2),v_top_zxy(plot_mask_top,3),v_top_zxy(plot_mask_top,1),'r.')
hold on
scatter3(v_btm_zxy(plot_mask_btm,2),v_btm_zxy(plot_mask_btm,3),v_btm_zxy(plot_mask_btm,1),'b.')
axis equal

stfig('Alignment of the BECs test');
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
clf
plot(top_halo.shot_num,...
    top_halo.num_counts,...
    'kx-','LineWidth',1.5)
hold on
plot(bottom_halo.shot_num,...
    bottom_halo.num_counts,...
    'bx-','LineWidth',1.5)
legend('top halo','bottom halo')
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

stfig('Vel Width History test');
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
subplot(2,3,1)
ndhist(v_btm_zxy(:,2:3));
hold on
plot(cos(t),sin(t),'k')
axis equal
xlabel('$v_x$')
ylabel('$v_y$')
subplot(2,3,2)
ndhist(v_btm_zxy(:,[1,3]));
hold on
plot(cos(t),sin(t),'k')
axis equal
xlabel('$v_z$')
ylabel('$v_y$')
subplot(2,3,3)
ndhist(v_btm_zxy(:,1:2));
hold on
plot(cos(t),sin(t),'k')
axis equal
xlabel('$v_z$')
ylabel('$v_x$')

subplot(2,3,4)
ndhist(v_btm_zxy_unnorm(:,2:3));
axis equal
xlabel('$v_x$')
ylabel('$v_y$')
subplot(2,3,5)
ndhist(v_btm_zxy_unnorm(:,[1,3]));
axis equal
xlabel('$v_z$')
ylabel('$v_y$')
subplot(2,3,6)
ndhist(v_btm_zxy_unnorm(:,1:2));
axis equal
xlabel('$v_z$')
ylabel('$v_x$')

stfig('average radius vs angle')
clf
nbins=50;
theta_bins = linspace(-pi,pi,nbins);
phi_bins = linspace(-pi/2,pi/2,nbins);
[theta_top,~] = cart2pol(v_top_zxy(:,2),v_top_zxy(:,3));
[theta_btm,~] = cart2pol(v_btm_zxy(:,2),v_btm_zxy(:,3));
phi_top = atan(v_top_zxy(:,1)./sqrt(v_top_zxy(:,2).^2+v_top_zxy(:,3).^2));
phi_btm = atan(v_btm_zxy(:,1)./sqrt(v_btm_zxy(:,2).^2+v_btm_zxy(:,3).^2));
for ii = 1:(nbins-1)
    r_btm_zxy_masked = r_dist_btm_unnorm(theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1));
    r_top_zxy_masked = r_dist_top_unnorm(theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1));
    v_btm_r(ii,1) = mean(r_btm_zxy_masked);
    v_top_r(ii,1) = mean(r_top_zxy_masked);
    
    r_btm_zxy_masked = r_dist_btm_unnorm(phi_bins(ii)<phi_btm & phi_btm<=phi_bins(ii+1));
    r_top_zxy_masked = r_dist_top_unnorm(phi_bins(ii)<phi_top & phi_top<=phi_bins(ii+1));
    v_btm_r(ii,2) = mean(r_btm_zxy_masked);
    v_top_r(ii,2) = mean(r_top_zxy_masked);
    
    theta(ii) = mean(theta_bins(ii:(ii+1)));
    phi(ii) = mean(phi_bins(ii:(ii+1)));
end
subplot(2,1,1)
plot(theta,v_btm_r(:,1),'linewidth',1.5)
hold on
plot(theta,v_top_r(:,1),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\theta\)')
subplot(2,1,2)
plot(phi,v_btm_r(:,2),'linewidth',1.5)
hold on
plot(phi,v_top_r(:,2),'linewidth',1.5)
legend('bottom','top')
ylabel('Average radial value')
xlabel('\(\phi\)')

stfig('spherical density plot')
for ii = 1:(nbins-1)
    for jj = 1:(nbins-1)
        ang_mask_btm = theta_bins(ii)<theta_btm & theta_btm<=theta_bins(ii+1) ...
            & phi_bins(jj)<phi_btm & phi_btm<=phi_bins(jj+1);
        area = abs(theta_bins(ii+1)-theta_bins(ii))*abs(sin(phi_bins(jj+1))-sin(phi_bins(jj)));
        ang_mask_top = theta_bins(ii)<theta_top & theta_top<=theta_bins(ii+1) ...
            & phi_bins(jj)<phi_top & phi_top<=phi_bins(jj+1);
        density_btm(ii,jj) = sum(ang_mask_btm)./area;
        avg_r_btm(ii,jj)  = nanmean(r_dist_btm_unnorm(ang_mask_btm));
        density_top(ii,jj) = sum(ang_mask_top)./area;
        avg_r_top(ii,jj)  = nanmean(r_dist_top_unnorm(ang_mask_top));
    end
end
subplot(2,2,1)
pcolor(phi,theta,avg_r_btm)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('average radius bottom')
subplot(2,2,2)
title('average density bottom')
pcolor(phi,theta,density_btm)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('density')
subplot(2,2,3)
pcolor(phi,theta,avg_r_top)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('average radius top')
subplot(2,2,4)
title('average density top')
pcolor(phi,theta,density_top)
shading flat
colorbar
ylabel('\(\theta\)')
xlabel('\(\phi\)')
title('density')

%%
stfig('radial density');
clf
rad_shift=0.059;
z_shift_top = [rad_shift.*ones(size(v_top_zxy_unnorm,1),1),zeros(size(v_top_zxy_unnorm,1),2)];
z_shift_btm = [rad_shift.*ones(size(v_btm_zxy_unnorm,1),1),zeros(size(v_btm_zxy_unnorm,1),2)];


[theta_top,rxy_top] = cart2pol(v_top_zxy_unnorm(:,2),v_top_zxy_unnorm(:,3));
combined_vzr = [v_top_zxy_unnorm(:,1)+z_shift_top(:,1),rxy_top;...
        v_top_zxy_unnorm(:,1)+z_shift_top(:,1),-rxy_top].*1e3;
    
ndhist(combined_vzr(:,[2,1]),'bins',4,'filter');
%caxis([0 4])
    hold on
    plot(zeros(1,1000),linspace(-0.14,0.14,1000).*1e3,'k-','LineWidth',3.8)
    xlabel('$v_r$ (mm/s)')
    %ylabel('$v_z$ (mm/s)')
    colormap('default')
    axis equal
    %caxis([0 9])
    ylim([-0.13,0.12].*1e3)
    ax = gca;
k = 0.02;%
ax.TickLength = [k, k]; % Make tick marks longer.
ax.LineWidth = 100*k; % Make tick marks thicker.

end