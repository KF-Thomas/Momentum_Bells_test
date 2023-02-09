function output = halo_reanalysis(halo_counts_vel, logging)
    arguments
        halo_counts_vel 
        logging = true
    end 
    % halo_reanalysis.m
    %   wraps the output from halo_analysis
    %   still in development (09 Feb 2023 - Tony)
    %     
    %{
        halo_reanalysis(halo{1}.counts_vel);
        halo_reanalysis(halo{1}.counts_vel, false);
    %}

    %%  Convert the dataset into a more usable format (for fitting purpose)

    shots_total = size(halo_counts_vel,1);
    %{
    xy_proj_normed = cell(shots_total, 1);
    for is = 1:shots_total
        this_shot_halo = halo_counts_vel{is};
        x = this_shot_halo(:, 2);
        y = this_shot_halo(:, 3);
        z = this_shot_halo(:, 1);
        
        r_xy = sqrt(x.^2 + y.^2 + z.^2);
        r_xy = 1;
        xn = x ./ r_xy;
        yn = y ./ r_xy;

        xy_proj_normed{is} = [xn yn];
        
    end 
    clear x y z ;
    xy_proj_normed_flat = cell2mat(xy_proj_normed);
    %}
    
%     figure(230);clf(230);figure(230);
%     scatter(xy_proj_normed_flat(:,1), xy_proj_normed_flat(:,2) ,'.');
%     view(0,90);
    
    halo_zxy = cell2mat(halo_counts_vel);
    halo_xyz = [halo_zxy(:,2), halo_zxy(:,3), halo_zxy(:,1)];

%     lim_r = [0.04, 0.08];
%     lim_


    %% plot xy projection
    figure(231);clf(231);figure(231);
    scatter3(halo_xyz(:,1), halo_xyz(:,2), halo_xyz(:,3) , ...
        '.','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1, 'LineWidth',0.1); hold on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title("halo projection to xy");
    view(0,90);
    xlim([-0.1, 0.1]);
    ylim([-0.1, 0.1]);
    zlim([-0.3, 0.3]);
    
    ori = [0 0 0];
    % https://au.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
%     [ecp, erp, evp, epp, ehp] = ellipsoid_fit(halo_xyz, '0');
    [ecg, erg, evg, epp, ehg] = ellipsoid_fit(halo_xyz);

    if logging 
        fprintf("    center,   radius     initial \n");
        disp([ecg, erg]);
    end 

    halo_xyz_c = [halo_xyz(:,1)-ecg(1), halo_xyz(:,2)-ecg(2), halo_xyz(:,3)-ecg(3)];
%     [ecg, erg, evg, epp, ehg] = ellipsoid_fit(halo_xyz_c);
    
    ievg = inv(evg);
    halo_xyz_i = (ievg * halo_xyz_c')' ; 
    halo_xyz_i_n = [halo_xyz_i(:,1)/erg(1), halo_xyz_i(:,2)/erg(2), halo_xyz_i(:,3)/erg(3)];
    halo_xyz_n = (evg * halo_xyz_i_n')';

    rotate3d(gcf, 'on');


    %% plot halo ellipsoid fit
    figure(233);clf(233);figure(233);
%     eev = evg * mean(erg);
%     plot3([ecg(1) eev(1,1)],[ecg(2) eev(1,2)],[ecg(3) eev(1,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(2,1)],[ecg(2) eev(2,2)],[ecg(3) eev(2,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(3,1)],[ecg(2) eev(3,2)],[ecg(3) eev(3,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(1,1)],[ecg(2) eev(1,2)],[ecg(3) eev(1,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(2,1)],[ecg(2) eev(2,2)],[ecg(3) eev(2,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(3,1)],[ecg(2) eev(3,2)],[ecg(3) eev(3,3)], 'r-^','LineWidth',3);

%     legend(('first', 'second', 'third'})
  
%     ellipsoid(ecp(1), ecp(2), ecp(3), erp(1), erp(2), erp(3),100); 
%     alpha 0.5;
%     grid off;
%     mesh off;
%     ep.EdgeColor = 'none';
%     set(ep, 'EdgeColor', 'none');

%     xyz_hacked = [(halo_xyz(:,1)-ecp(1))/erp(1),(halo_xyz(:,2)-ecp(2))/erp(2),(halo_xyz(:,3)-ecp(3))/erp(3)];
%     
%     [ech, erh, evh, eph, ehh] = ellipsoid_fit(xyz_hacked, '0');
% %     disp({ech, erh, evh, eph, ehh});
% 
%     figure(232);clf(232);figure(232);
%     scatter3(xyz_hacked(:,1), xyz_hacked(:,2), xyz_hacked(:,3) ,'.'); hold on;
%     ellipsoid(ech(1), ech(2), ech(3), erh(1), erh(2), erh(3),100); 

    
    [ecc, erc, evc, epc, ehc] = ellipsoid_fit(halo_xyz_n,'0');
    ellipsoid(ecc(1), ecc(2), ecc(3), erc(1), erc(2), erc(3),100);  hold on;
%     alpha 0.5 ;
%     alpha(sc, 1); % I hate matlab... this is so easy to do in plt.plot 
    findobj_result = findobj(gcf,'type','Surface'); % warning: I'm hard coding this plot

    alpha(findobj_result, 0.4); % WTF why doesn't this change the black lines?! I hate matlab 
    set(findobj_result, 'EdgeColor','none'); % WTF I can't set color, this can only be on/off??? I hate matlab

    scatter3(halo_xyz_n(:,1), halo_xyz_n(:,2), halo_xyz_n(:,3) ,'.'); hold on;

    xlabel("x");
    ylabel("y");
    zlabel("z");
    title("Fitting halo to ellipsoid")
    view(0,90);
    rotate3d(gcf, 'on');

    xlim([-1,1]*1.2);
    ylim([-1,1]*1.2);
    zlim([-1,1]*1.2);

    if logging
        fprintf("    center,   radius     final \n");
        disp([ecc, erc]);
    end 
    hold off;
    

    %%  distribution around halo azm 
    %   check uniformity 

%     azm_hist_bin = 360;
    [halo_radial,halo_azm,halo_elev]=ConvToSph(halo_zxy);
    figure(234);clf(234);figure(234);
    
    histogram(halo_azm,30 ,'Normalization','pdf','DisplayName',' 30 bins', ...
        'FaceColor','blue','FaceAlpha',0.3,'EdgeColor','none'); hold on; 
    histogram(halo_azm,360,'Normalization','pdf','DisplayName','360 bins', ...
        'FaceColor','red', 'FaceAlpha',0.3,'EdgeColor','none'); hold on;
    hist_counts = histcounts(halo_azm, 360, 'Normalization','pdf');
    
    legend('AutoUpdate','off');
    xlabel("azimuthal");
    ylabel("counts");

    halo_azm_mean = mean(hist_counts);
    halo_azm_std  = std(hist_counts);
    
    rectangle('Position',[0,halo_azm_mean-halo_azm_std 2*pi 2*halo_azm_std], ...
        'LineStyle','none','FaceColor',[0, 1, 0,.2]); hold on;
    yline(halo_azm_mean, '-',num2str(halo_azm_mean,4)+"$\pm$"+num2str(halo_azm_std,4), ...
        'LineWidth',3,'Color',[0.1, 0.8, 0.1, 1],'Interpreter','latex'); hold on; 
    
    yline(1/(2*pi), '-','$\frac{1}{2\pi}$','LabelVerticalAlignment','bottom', ...
        'LineWidth',1,'Color','cyan','Interpreter','latex'); hold on; 
    
    
%%%%    235 
    figure(235);clf(235);figure(235);
    
    histogram2(halo_azm, halo_elev, [90,45], 'FaceColor','flat','ShowEmptyBins','on', 'FaceAlpha',0.95,'EdgeColor','none'); hold on;

    % set(gca, 'Projection','perspective');
    xlabel("azimuthal");
    ylabel("elevation");
    zlabel("counts");
    
    view(-5,85);
    rotate3d(gcf, 'on');
    %set(gcf, "GraphicsSmoothing", "on");


%%%%    236
    figure(236);clf(236);figure(236);
    
    halo_azm_0pi = halo_azm(halo_azm<pi);
    halo_azm_pi2 = halo_azm(halo_azm>=pi) - pi ;

    histogram(halo_azm_pi2,180,'Normalization','pdf','DisplayName','0 to $\pi$', ...
        'FaceColor','red', 'FaceAlpha',0.3,'EdgeColor','none'); hold on;
    histogram(halo_azm_0pi,180,'Normalization','pdf','DisplayName','$\pi$ to $2\pi$', ...
        'FaceColor','blue', 'FaceAlpha',0.3,'EdgeColor','none'); hold on;

    xlabel("azimuthal");
    ylabel("counts");
    
    legend();

%%%%    237
    figure(237);clf(237);figure(237);

    scatter(halo_azm, halo_radial, '.'); hold on;

    xlabel("azimuthal");
    ylabel("radial");
%     
%     alin = linspace(0,2*pi,2);
%     yL_lin   = -0.055*(alin-0.50)+0.058;
%     yL_lin_L = -0.055*(alin-0.50+0.1)+0.058;
%     yL_lin_R = -0.055*(alin-0.50-0.1)+0.058;
%     yR_lin   = +0.055*(alin-3.35)+0.058;
%     yR_lin_L = +0.055*(alin-3.35+0.1)+0.058;
%     yR_lin_R = +0.055*(alin-3.35-0.1)+0.058;
% 
%     plot(alin, yL_lin); hold on;
%     plot(alin, yR_lin); hold on;
%     plot(alin, yL_lin_L);
%     plot(alin, yL_lin_R);
%     plot(alin, yR_lin_L);
%     plot(alin, yR_lin_R);
% 
    xlim([0 2]*pi);
    ylim([0.05,0.08]);
%     hold off;

%%%%    238
    figure(238);clf(238);figure(238);

    scatter3(halo_azm, halo_elev, halo_radial, '.');

    xlabel("azimuthal");
    ylabel("elevation");
    zlabel("radial")
    xlim([0 2]*pi);
    ylim([0 1]*pi);
%     zlim([0.05 0.1]);

    view(0,90);
    rotate3d(gcf, 'on');

%%%%    239
    figure(239);clf(239);figure(239);

    scatter3(halo_azm, halo_elev, halo_radial, '.');

    xlabel("azimuthal");
    ylabel("elevation");
    zlabel("radial")
    xlim([0 2]*pi);
    ylim([0 1]*pi);
    zlim([0.05 0.08]);

    view(90,0);
    rotate3d(gcf, 'on');


    shg

    %% for the sake of cross compatability
%     output = cellfun(@(shot) [(shot(:,1)-ecp(1))/erp(1),(shot(:,2)-ecp(2))/erp(2),(shot(:,3)-ecp(3))/erp(3)] ,halo_counts_vel, 'UniformOutput',false);

    output = cell(shots_total,1);
    for is = 1:shots_total
        this_shot_halo = halo_counts_vel{is};
        x = this_shot_halo(:, 2);
        y = this_shot_halo(:, 3);
        z = this_shot_halo(:, 1);
    
        [hr,ha,he] = ConvToSph([z, x, y]);
        
        cutL = (hr > -0.055*(ha-0.50+0.1)+0.058) & (hr < -0.055*(ha-0.50-0.1)+0.058);
        cutR = (hr < +0.055*(ha-3.35+0.1)+0.058) & (hr > +0.055*(ha-3.35-0.1)+0.058);
        cut = cutL | cutR; 

        xyz = [x, y, z];
%         xyz = xyz(~cut, :);

        xyz_i = (ievg * xyz')' ; 
        xyz_i_n = [xyz_i(:,1)/erg(1), xyz_i(:,2)/erg(2), xyz_i(:,3)/erg(3)];
        xyz_n = (evg * xyz_i_n')';
        
        zxy = [xyz_n(:,3), xyz_n(:,1), xyz_n(:,2)];
        
        output{is} = zxy;
    end 


end










