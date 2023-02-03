function output = halo_reanalysis(halo_counts_vel)
    arguments
        halo_counts_vel 

    end 

    % halo_re_out = halo_reanalysis(halo{1}.counts_vel); 

    shots_total = size(halo_counts_vel,1);

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
    
    
%     figure(230);clf(230);figure(230);
%     scatter(xy_proj_normed_flat(:,1), xy_proj_normed_flat(:,2) ,'.');
    
    
    halo_counts_vel_flat = cell2mat(halo_counts_vel);
    halo_xyz = [halo_counts_vel_flat(:,2), halo_counts_vel_flat(:,3), halo_counts_vel_flat(:,1)];
%     view(0,90);

    figure(231);clf(231);figure(231);
    scatter3(halo_xyz(:,1), halo_xyz(:,2), halo_xyz(:,3) ,'.'); hold on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    view(0,90);
    

    ori = [0 0 0];

    % https://au.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
%     [ecp, erp, evp, epp, ehp] = ellipsoid_fit(halo_xyz, '0');
    [ecg, erg, evg, epp, ehg] = ellipsoid_fit(halo_xyz);
    disp("initial center, radius")
    disp([ecg, erg]);

    halo_xyz_c = [halo_xyz(:,1)-ecg(1), halo_xyz(:,2)-ecg(2), halo_xyz(:,3)-ecg(3)];
    [ecg, erg, evg, epp, ehg] = ellipsoid_fit(halo_xyz_c);
    
    ievg = inv(evg);
    halo_xyz_i = (ievg * halo_xyz_c')' ; 
    halo_xyz_i_n = [halo_xyz_i(:,1)/erg(1), halo_xyz_i(:,2)/erg(2), halo_xyz_i(:,3)/erg(3)];
    halo_xyz_n = (evg * halo_xyz_i_n')';

    figure(233);clf(233);figure(233);
    scatter3(halo_xyz_n(:,1), halo_xyz_n(:,2), halo_xyz_n(:,3) ,'.'); hold on;


%     eev = evg * 0.07;
%     plot3([ecp(1) eev(1,1)],[ecp(2) eev(1,2)],[ecp(3) eev(1,3)], 'r-^','LineWidth',3);
%     plot3([ecp(1) eev(2,1)],[ecp(2) eev(2,2)],[ecp(3) eev(2,3)], 'r-^','LineWidth',3);
%     plot3([ecp(1) eev(3,1)],[ecp(2) eev(3,2)],[ecp(3) eev(3,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(1,1)],[ecg(2) eev(1,2)],[ecg(3) eev(1,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(2,1)],[ecg(2) eev(2,2)],[ecg(3) eev(2,3)], 'r-^','LineWidth',3);
%     plot3([ecg(1) eev(3,1)],[ecg(2) eev(3,2)],[ecg(3) eev(3,3)], 'r-^','LineWidth',3);
    
    xlabel("x");
    ylabel("y");
    zlabel("z");

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
    ellipsoid(ecc(1), ecc(2), ecc(3), erc(1), erc(2), erc(3),100); 
    alpha 0.1;

    xlabel("x");
    ylabel("y");
    zlabel("z");

    view(0,90);
    disp("final center, radius")
    disp([ecc, erc]);

%     output = cellfun(@(shot) [(shot(:,1)-ecp(1))/erp(1),(shot(:,2)-ecp(2))/erp(2),(shot(:,3)-ecp(3))/erp(3)] ,halo_counts_vel, 'UniformOutput',false);

    output = cell(shots_total,1);
    for is = 1:shots_total
        this_shot_halo = halo_counts_vel{is};
        x = this_shot_halo(:, 2);
        y = this_shot_halo(:, 3);
        z = this_shot_halo(:, 1);
        xyz = [x, y, z];
        xyz_i = (ievg * xyz')' ; 
        xyz_i_n = [xyz_i(:,1)/erg(1), xyz_i(:,2)/erg(2), xyz_i(:,3)/erg(3)];
        xyz_n = (evg * xyz_i_n')';
        
        zxy = [xyz_n(:,3), xyz_n(:,1), xyz_n(:,2)];
        
        output{is} = zxy;
    end 


end







