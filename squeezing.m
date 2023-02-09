function out=squeezing(halo_centered_cells,isverbose)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%new file 1/2/23 - M.W

%calulates squeezing
%inputs 
    %cell array for eahc file containing Z,X,Y counts
    %isverbose = 1 displays some info, =0 does not. 
%outputs
    %a cell that contains two matrices, the first contains what is needed
    %for the sqz angle plot (angle,variance,variance_se)
    %the second which contains [mean_opst_bin,unc_opst_bin,mean_other_bin,unc_other_bin]
%this takes the centered halo data as cells for each shot of z,x,y 
%converts it to radial cord and then finds the normalized number vairance
%the traditional way (http://arxiv.org/pdf/1008.0845.pdf)of 
%showing this squezing data is to plot the disjointed pairs and the jointed
%pairs, this seems crude and arbitrary( the x axis has no information)
%instead i plot using the angle between the bin centers when this is pi
%then you expect decreased num var above it not

% ref system
% BEC's are at the poles
% azimuthal angle is arround the equator (zero to 2 pi)
% inclination is angle from the poles (which is zero inclination)(zero to pi)
% this was chosen to reduce the wrapping problems arround zero/2 pi 

%TO FIX
%cuts in elev direction are currently unsupported 9/2/23

%To add
%artificial QE to look at how QE affects squeezing

%================================ START USER INPUT=======================================

%----------- Plots/Calculations
plot_sqz_bins=0; %plot sqz with bin #
plot_sqz_angle=1; %plot sqz with angle

plot_sqz_num_zones=1; %plot average norm var as number of bins changes. Returns QE estimate from fit.

plot_sqz_counts=1; %plot sqz with angle for shots with different number of total counts.
%Data above, below and within range count_lims are plotted as separate categories.
%window_counts must be set to zero for this //not anymore bitches

%----------- Data 
window_counts=1; %Only use shots with counts in range count_lims
count_lims=[15,60]; %acceptable range for counts
plot_count_distribution=1; %Plot histogram of total counts in each shot

%Regions to cut out of halo (e.g. if there are artefacts)
%To implement a cut between azm angles az1,az2 and elev angles elev1,elev2:
% cuts_azm = [az1 az2; ...]
% cuts_elev = [elev1,elev2; ...]
cuts_azm= [];%[1.3 1.45]*pi; % [];
cuts_elev= [];%[0 1]*pi;   % []; %currently doesn't support elev cuts that don't go entirely from 0 to pi 
remove_opposite_zones=1; %for each cut region, remove the opposite zone

%----------- Define Bins
mirror_azm=0; %use for the mag sensitive halo which you must cut out a region

range_azm=[0 2]*pi;%azimuthal range for bins. to look at whole halo, [0 2]*pi
range_elev= [0 1]*pi; %[0.0 1]*pi;%Elevation range. whole halo [0.0 1]*pi

steps_azm=50; %number of bins in azimuthal direction
steps_elev=2; %number of bins in elev direction

%defining the width seprately allows for under or over sampling bins
bin_width_azm=2*range(range_azm)/steps_azm; 
bin_width_elev=range(range_elev)/steps_elev;

%----------- Plot Halo Options
%halo plots require an even number of azm bins
plot_halo_zones=1; %plot halo with one corresponding zone pair highlighted
plot_all_halo_zones=1; %plot halo with all zones. If zones overlap, this will only show up for the last zone pair to be plotted.

plot_cut_lines=1; %plot cut lines on halo_zones plot
plot_with_cut=1;  %remove cut data from plot
plot_all_with_cut=1; %remove cut data from all_halo_zones plot
plot_all_with_cut_lines=0; %plot cut lines on halo_zones plot

%===================================END USER INPUT========================================

bin_info=struct('range_azm',range_azm,'range_elev',range_elev,...
        'steps_azm',steps_azm,'steps_elev',steps_elev,...
        'bin_width_azm',bin_width_azm,'bin_width_elev',bin_width_elev,...
        'cuts_azm',cuts_azm,'cuts_elev',cuts_elev,...
        'mirror_azm',mirror_azm);
output_bins = squeezing_bins(halo_centered_cells,bin_info,0);
output_norm_var = squeezing_norm_var(output_bins{1},output_bins{2},output_bins{3},window_counts,count_lims,plot_count_distribution,0);

angle_var_sd = output_norm_var{1};
mean_info = output_norm_var{2};

mean_opst_bin=mean_info(1,1);
unc_opst_bin=mean_info(1,2);
min_opst_bin=mean_info(1,3);
mean_other_bin=mean_info(2,1);
unc_other_bin=mean_info(2,2);
min_other_bin=mean_info(2,3);

if  isverbose
    disp(['mean opst bin ' ,num2str(mean_opst_bin),'±',num2str(unc_opst_bin)])
    disp(['mean other bins ',num2str(mean_other_bin),'±',num2str(unc_other_bin)])
    disp(['min opst bin ',num2str(min_opst_bin)])
end

%--------------------------------------------------------------------------

if plot_sqz_angle
    
    figure(10)
    errorbar(angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3),'x')
    %title(['Squeezing (',num2str(steps_azm),' azm ',num2str(steps_elev),' elev )']);
    set(gcf,'Color',[1 1 1]);
    xlabel('Angle Between Bins /Pi')
    ylabel('Normalised Variance')
    line([-2 2], [1 1],'Color','red');
    set(gca,'Xlim',[-0.1 1.1])

end

%--------------------------------------------------------------------------

if plot_sqz_bins

    norm_var_pair=output_norm_var{3};
    norm_var_opst=norm_var_pair{1};
    norm_var_rest=norm_var_pair{2};
    
    figure(12)
    plot(norm_var_rest(:,1),norm_var_rest(:,2),'+',norm_var_opst(:,1),norm_var_opst(:,2), 'x')
    %plot(norm_var_pair_num_rest(:,1),norm_var_pair_num_rest(:,2),'+',norm_var_pair_num_opst(:,1),norm_var_pair_num_opst(:,2), 'x')
    xlabel('Bin Pair Number')
    ylabel('Normalised Variance')
    line([0 size(norm_var_rest,1)], [1 1],'Color','red');
end

%--------------------------------------------------------------------------

if plot_sqz_num_zones    
    %define array with number of zones Nz
    bin_num_array = 10:10:400;
    
    if isverbose
        disp(['calculating for different numbers of zones'])
    end 
    
    %calculate mean normalised variance for correlated/uncorrelated zones
    %for each Nz
    %this could be a lot faster if we sort into bins first and then
    %combine bins, instead of rebinning data every time
    mean_var_output = [];
    bin_info_Nz=struct('range_azm',range_azm,'range_elev',range_elev,...
            'steps_azm',steps_azm,'steps_elev',steps_elev,...
            'bin_width_azm',bin_width_azm,'bin_width_elev',bin_width_elev,...
            'cuts_azm',cuts_azm,'cuts_elev',cuts_elev,...
            'mirror_azm',mirror_azm);
        
    for i=1:length(bin_num_array)
        n=bin_num_array(i);      

        bin_info_Nz.steps_azm=n;
        bin_info_Nz.bin_width_azm=(2*range(range_azm)/n);
        
        output_bins2 = squeezing_bins(halo_centered_cells,bin_info_Nz,0);
        output = squeezing_norm_var(output_bins2{1},output_bins2{2},output_bins2{3},window_counts,count_lims,0,0);
        %output = squeezing_bins(halo_centered_cells,range_azm, range_elev,n, steps_elev, (2*range(range_azm)/n), bin_width_elev,cuts_azm,cuts_elev,1,count_lims,0);
        output2 = output{2};
        mean_var_output = [mean_var_output; output2(1,:) output2(2,:)];
    end
    
    %fit 
    %select fit region
    fit_region_start=1;
    fit_region_end=length(bin_num_array);
    
    %function to fit: V=1-QE*P(Nz)
    halo_radius=0.065;
    halo_circumference=halo_radius*2*pi;
    var_function = @(arg, nn) 1 - arg(1).*(nn./halo_circumference).*((-1+exp(-halo_circumference^2./(2*(nn.*arg(2)).^2))).*sqrt(2/pi).*arg(2) + ...
               (halo_circumference./nn).*erf(halo_circumference./(sqrt(2).*nn.*arg(2))));
    beta0 = [0.1, halo_radius*0.05];
    mdl = fitnlm(bin_num_array(fit_region_start:fit_region_end),mean_var_output(fit_region_start:fit_region_end,1),var_function,beta0);
    if isverbose
        disp(['QE estimate from fit: ',num2str(mdl.Coefficients.Estimate(1)),'±',num2str(mdl.Coefficients.SE(1))])
        disp(['$\sigma$ (~position uncertainty): ',num2str(mdl.Coefficients.Estimate(2)),'±',num2str(mdl.Coefficients.SE(2))])
    end
    
    %curve of best fit to plot
    nz_lin=linspace(bin_num_array(1),bin_num_array(end),1000)';
    [y_mod, y_mod_ci] = predict(mdl, nz_lin);
    
    %Plot
    stfig('Variance and number of zones');
    clf
    errorbar(bin_num_array,mean_var_output(:,1),mean_var_output(:,2),'.','DisplayName',"Correlated Zones")
    hold on     
    errorbar(bin_num_array,mean_var_output(:,4),mean_var_output(:,5),'.','DisplayName',"Uncorrelated Zones")
    
    plot(nz_lin, y_mod,"Color",'r', "DisplayName", "Fit");
    fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], 'r', ...
        'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Fit $2\sigma$ CI")
    
    line([0 bin_num_array(end)], [1 1],'Color','red','DisplayName','1');
    xlabel('Number of Zones')
    ylabel('Normalised Variance')
    legend();
    
end

%--------------------------------------------------------------------------

if plot_sqz_counts==1
    
    bin_pairs=output_bins{1};
    total_counts=output_bins{2};
    
    count_lims_ok = count_lims;
    count_lims_low = [0,count_lims(1)];
    count_lims_high = [count_lims(2),inf];
    
    output_norm_var_ok = squeezing_norm_var(output_bins{1},output_bins{2},output_bins{3},1,count_lims_ok,0,0);
    angle_var_sd_ok=output_norm_var_ok{1};
    output_norm_var_low = squeezing_norm_var(output_bins{1},output_bins{2},output_bins{3},1,count_lims_low,0,0);
    angle_var_sd_low=output_norm_var_low{1};
    output_norm_var_high = squeezing_norm_var(output_bins{1},output_bins{2},output_bins{3},1,count_lims_high,0,0);
    angle_var_sd_high=output_norm_var_high{1};
    
    %plot stuff
    stfig('Norm Var vs Angle between bins (number of total counts/shot)');
    clf
    errorbar(angle_var_sd_ok(:,1),angle_var_sd_ok(:,2),angle_var_sd_ok(:,3),'x')
    hold on 
    errorbar(angle_var_sd_low(:,1),angle_var_sd_low(:,2),angle_var_sd_low(:,3),'x')
    errorbar(angle_var_sd_high(:,1),angle_var_sd_high(:,2),angle_var_sd_high(:,3),'x')
    hold off
    %title(['Squeezing (',num2str(steps_azm),' azm ',num2str(steps_elev),' elev )']);
    set(gcf,'Color',[1 1 1]);
    xlabel('Angle Between Bins /Pi')
    ylabel('Normalised Variance')
    line([-2 2], [1 1],'Color','red');
    set(gca,'Xlim',[-0.1 1.1])
    legend('ok','low','high')
end 

%--------------------------------------------------------------------------

if plot_halo_zones
    
    %halo data to plot
    v_zxy = cell2mat(halo_centered_cells');
    plot_mask = rand(size(v_zxy,1),1)<0.65; %controls how sparse the plot is
    v_x=v_zxy(plot_mask,2);
    v_y=v_zxy(plot_mask,3);
    v_z=v_zxy(plot_mask,1);
    
    %create bins
    bins = create_bins(bin_info);
    bin_centers_azm = bins{1};
    bin_pairs_azm = bins{2};
    bin_centers_elev = bins{3};
    bin_pairs_elev = bins{4}; 
   
    %pick out two corresponding zones (one given by m,p)
    m=1 ;
    p=1 ;
    m2 = mod(m+steps_azm/2, steps_azm);
    if m2==0
        m2 = steps_azm;
    end
    p2=steps_elev-p+1;
    zone_1_azm = bin_pairs_azm(m,:);
    zone_1_elev = bin_pairs_elev(p,:);
    zone_2_azm = bin_pairs_azm(m2,:); %assume steps_azm even
    zone_2_elev = bin_pairs_elev(p2,:);
    
    %convert to spherical coordinates - matches ConvToSph.m
    r = sqrt(v_x.^2+v_y.^2+v_z.^2);
    halo_azm = atan2(v_y,v_x)+pi;
    halo_elev = acos(v_z./r);
    
    %cuts
    if plot_with_cut==1
        for i=1:size(cuts_azm,1)
            cut = halo_azm>cuts_azm(i,1) & halo_azm<cuts_azm(i,2) & halo_elev>cuts_elev(i,1) & halo_elev<cuts_elev(i,2);
            if remove_opposite_zones==1
             cut2 = halo_azm>wrapTo2Pi(cuts_azm(i,1)+pi) & halo_azm<wrapTo2Pi(cuts_azm(i,2)+pi) & halo_elev<(pi-cuts_elev(i,1)) & halo_elev>(pi-cuts_elev(i,2));
             cut = or(cut, cut2);
            end
            r=r(~cut);
            halo_azm=halo_azm(~cut);
            halo_elev = halo_elev(~cut);
            v_x=v_x(~cut);
            v_y=v_y(~cut);
            v_z=v_z(~cut);
        end
    end
    
    %mask for each zone
    %wrap1 and wrap2 masks account for bins wrapping over 0,2pi
    zone_1_mask = zone_1_azm(1)<zone_1_azm(2)...
        & halo_azm>zone_1_azm(1) & halo_azm<zone_1_azm(2) & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);
    zone_1_mask_wrap1 = zone_1_azm(1)>zone_1_azm(2)...
        & halo_azm>zone_1_azm(1) & halo_azm<2*pi & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);
    zone_1_mask_wrap2 = zone_1_azm(1)>zone_1_azm(2)...
        & halo_azm>0 & halo_azm<zone_1_azm(2) & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);
    
    zone_2_mask = zone_2_azm(1)<zone_2_azm(2)...
        & halo_azm>zone_2_azm(1) & halo_azm<zone_2_azm(2) & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);
    zone_2_mask_wrap1 = zone_2_azm(1)>zone_2_azm(2)...
        & halo_azm>zone_2_azm(1) & halo_azm<2*pi & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);
    zone_2_mask_wrap2 = zone_2_azm(1)>zone_2_azm(2)...
        & halo_azm>0 & halo_azm<zone_2_azm(2) & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);
    
    %halo data for each zone
    halo_zone_1_x = [v_x(zone_1_mask); v_x(zone_1_mask_wrap1); v_x(zone_1_mask_wrap2)];
    halo_zone_1_y = [v_y(zone_1_mask); v_y(zone_1_mask_wrap1); v_y(zone_1_mask_wrap2)];
    halo_zone_1_z = [v_z(zone_1_mask); v_z(zone_1_mask_wrap1); v_z(zone_1_mask_wrap2)];
    
    halo_zone_2_x = [v_x(zone_2_mask); v_x(zone_2_mask_wrap1); v_x(zone_2_mask_wrap2)];
    halo_zone_2_y = [v_y(zone_2_mask); v_y(zone_2_mask_wrap1); v_y(zone_2_mask_wrap2)];
    halo_zone_2_z = [v_z(zone_2_mask); v_z(zone_2_mask_wrap1); v_z(zone_2_mask_wrap2)];
    
    %plot halo and zones
    stfig(['halo zones: (az,elev)= (',num2str(m),',',num2str(p),') and (',num2str(m2),',',num2str(p2),')']);
    clf
    
    scatter3(v_x,v_y,v_z,'.')
    hold on
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    zlabel('$v_z$')
    scatter3(halo_zone_1_x,halo_zone_1_y,halo_zone_1_z,'.','k')
    scatter3(halo_zone_2_x,halo_zone_2_y,halo_zone_2_z,'.','k')
    if plot_cut_lines==1 
        r_max = max(r);
        for i=1:size(cuts_azm,1)
            [x,y,z] = sph2cart(cuts_azm(i,1),cuts_elev(i,1),r_max);
            [x2,y2,z2] = sph2cart(cuts_azm(i,1),cuts_elev(i,2),r_max);
            [x3,y3,z3] = sph2cart(cuts_azm(i,2),cuts_elev(i,1),r_max);
            [x4,y4,z4] = sph2cart(cuts_azm(i,2),cuts_elev(i,2),r_max);
            line([x,x2], [y,y2], [z,z2], 'LineWidth', 1, 'Color', 'r');
            line([x3,x4], [y3,y4], [z3,z4], 'LineWidth', 1, 'Color', 'r');
        end
    end
    hold off    

end

%--------------------------------------------------------------------------

if plot_all_halo_zones
    
    %halo data to plot
    v_zxy = cell2mat(halo_centered_cells');
    plot_mask = rand(size(v_zxy,1),1)<0.65; %controls how sparse the plot is
    v_x=v_zxy(plot_mask,2);
    v_y=v_zxy(plot_mask,3);
    v_z=v_zxy(plot_mask,1);
    
    %convert to spherical coordinates - matches ConvToSph.m
    r = sqrt(v_x.^2+v_y.^2+v_z.^2);
    halo_azm = atan2(v_y,v_x)+pi;
    halo_elev = acos(v_z./r);
    
    %cuts
    if plot_all_with_cut==1
        %create bins
        bins = create_bins(bin_info);
        bin_centers_azm = bins{1};
        bin_pairs_azm = bins{2};
        bin_centers_elev = bins{3};
        bin_pairs_elev = bins{4};
        
        for i=1:size(cuts_azm,1)
            cut = halo_azm>cuts_azm(i,1) & halo_azm<cuts_azm(i,2) & halo_elev>cuts_elev(i,1) & halo_elev<cuts_elev(i,2);
            if remove_opposite_zones==1
             cut2 = halo_azm>wrapTo2Pi(cuts_azm(i,1)+pi) & halo_azm<wrapTo2Pi(cuts_azm(i,2)+pi) & halo_elev<(pi-cuts_elev(i,1)) & halo_elev>(pi-cuts_elev(i,2));
             cut = or(cut, cut2);
            end
            r=r(~cut);
            halo_azm=halo_azm(~cut);
            halo_elev = halo_elev(~cut);
            v_x=v_x(~cut);
            v_y=v_y(~cut);
            v_z=v_z(~cut);
        end
    end
    
    %plot    
    stfig('halo zones');
    clf
    
    hold on
    axis equal
    xlabel('$v_x$')
    ylabel('$v_y$')
    zlabel('$v_z$')
    
    zone_colours = lines(steps_azm*steps_elev/2);   
    for m=1:steps_azm/2
        for p=1:steps_elev
    
            %pick out two corresponding zones (one given by m,p)
            zone_1_azm = bin_pairs_azm(m,:);
            zone_1_elev = bin_pairs_elev(p,:);
            zone_2_azm = bin_pairs_azm(m+steps_azm/2,:); %assume steps_azm even %FIX THIS BIT
            zone_2_elev = bin_pairs_elev(steps_elev-p+1,:); %assume steps_elev even 

            %mask for each zone
            %wrap1 and wrap2 masks account for bins wrapping over 0,2pi
            zone_1_mask = zone_1_azm(1)<zone_1_azm(2)...
                & halo_azm>zone_1_azm(1) & halo_azm<zone_1_azm(2) & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);
            zone_1_mask_wrap1 = zone_1_azm(1)>zone_1_azm(2)...
                & halo_azm>zone_1_azm(1) & halo_azm<2*pi & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);
            zone_1_mask_wrap2 = zone_1_azm(1)>zone_1_azm(2)...
                & halo_azm>0 & halo_azm<zone_1_azm(2) & halo_elev>zone_1_elev(1) & halo_elev<zone_1_elev(2);

            zone_2_mask = zone_2_azm(1)<zone_2_azm(2)...
                & halo_azm>zone_2_azm(1) & halo_azm<zone_2_azm(2) & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);
            zone_2_mask_wrap1 = zone_2_azm(1)>zone_2_azm(2)...
                & halo_azm>zone_2_azm(1) & halo_azm<2*pi & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);
            zone_2_mask_wrap2 = zone_2_azm(1)>zone_2_azm(2)...
                & halo_azm>0 & halo_azm<zone_2_azm(2) & halo_elev>zone_2_elev(1) & halo_elev<zone_2_elev(2);

            %halo data for each zone
            halo_zone_1_x = [v_x(zone_1_mask); v_x(zone_1_mask_wrap1); v_x(zone_1_mask_wrap2)];
            halo_zone_1_y = [v_y(zone_1_mask); v_y(zone_1_mask_wrap1); v_y(zone_1_mask_wrap2)];
            halo_zone_1_z = [v_z(zone_1_mask); v_z(zone_1_mask_wrap1); v_z(zone_1_mask_wrap2)];

            halo_zone_2_x = [v_x(zone_2_mask); v_x(zone_2_mask_wrap1); v_x(zone_2_mask_wrap2)];
            halo_zone_2_y = [v_y(zone_2_mask); v_y(zone_2_mask_wrap1); v_y(zone_2_mask_wrap2)];
            halo_zone_2_z = [v_z(zone_2_mask); v_z(zone_2_mask_wrap1); v_z(zone_2_mask_wrap2)];
                        
            scatter3([halo_zone_1_x; halo_zone_2_x],[halo_zone_1_y; halo_zone_2_y],[halo_zone_1_z; halo_zone_2_z],'.','MarkerEdgeColor',zone_colours(m*p,:))
        end
    end
    
    if plot_all_with_cut_lines==1 
        r_max = max(r);
        for i=1:size(cuts_azm,1)
            [x,y,z] = sph2cart(cuts_azm(i,1),cuts_elev(i,1),r_max);
            [x2,y2,z2] = sph2cart(cuts_azm(i,1),cuts_elev(i,2),r_max);
            [x3,y3,z3] = sph2cart(cuts_azm(i,2),cuts_elev(i,1),r_max);
            [x4,y4,z4] = sph2cart(cuts_azm(i,2),cuts_elev(i,2),r_max);
            line([x,x2], [y,y2], [z,z2], 'LineWidth', 1, 'Color', 'r');
            line([x3,x4], [y3,y4], [z3,z4], 'LineWidth', 1, 'Color', 'r');
        end
    end
    
    hold off

end

out={[angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3)],[[mean_opst_bin,unc_opst_bin,min_opst_bin];[mean_other_bin,unc_other_bin,min_other_bin]]};

end

%delete

% %remove empty bins
% halo_centered_cells=halo_centered_cells(~cellfun('isempty',halo_centered_cells));

% % create the bins
% bin_centers_azm=wrapTo2Pi(linspace(range_azm(1),range_azm(2),steps_azm+1)-range(range_azm)/(2*steps_azm));
% bin_centers_azm=bin_centers_azm(2:end); %first and last elements are separated by 2pi and we only need one of them
% 
% %pairs containing start and end of each bin
% bin_pairs_azm=transpose([wrapTo2Pi(bin_centers_azm-bin_width_azm/2) ; wrapTo2Pi(bin_centers_azm+bin_width_azm/2)]);
% 
% 
% if mirror_azm % ? not sure what this does ?
%     bin_centers_azm=[bin_centers_azm, bin_centers_azm+pi];
%     bin_pairs_azm=[bin_pairs_azm;bin_pairs_azm+pi];
% end
% 
% bin_centers_elev=wrapTo2Pi(linspace(range_elev(1),range_elev(2),steps_elev+1)-range(range_elev)/(2*steps_elev));
% bin_centers_elev=bin_centers_elev(2:end);
% bin_pairs_elev=transpose([wrapTo2Pi(bin_centers_elev-bin_width_elev/2) ; wrapTo2Pi(bin_centers_elev+bin_width_elev/2)]);

% if isverbose
% disp('binning') %indicate that binning is complete
% end
% 
% % convert halo to spherical coordinates and count atoms in each bin
% % halo_counts=[]; %not used
% % halo_centered_spherical={}; %not used
% 
% shot_num=size(halo_centered_cells,2) %number of shots in halo
% if plot_sqz_counts||window_counts
%     total_counts=zeros(shot_num,1); %list to store total counts (indicating mode occupancy) of each halo
% end
% 
% angle_counts={}; %stores number of counts in each bin for each shot. Each element is an array containing the shot number, azm/elev bin centres and counts
% 
% for n=1:size(halo_centered_cells,2) %n=1:shot_num
%     angle_counts_file={};
%     single_halo=halo_centered_cells{n};
%     
%     if plot_sqz_counts||window_counts
%         total_counts(n)=size(single_halo,1);
%     end 
%     
%     [halo_radial,halo_azm,halo_elev]=ConvToSph(single_halo);
%     
%     %Remove data according to specified cuts
%     %don't need this anymore since bins have been cut
%     cut_data=0;
%     if cut_data
%      for i=1:size(cuts_azm,1)
%         cut = halo_azm>cuts_azm(i,1) & halo_azm<cuts_azm(i,2) & halo_elev>cuts_elev(i,1) & halo_elev<cuts_elev(i,2);
%         if remove_opposite_zones==1
%              cut2 = halo_azm>wrapTo2Pi(cuts_azm(i,1)+pi) & halo_azm<wrapTo2Pi(cuts_azm(i,2)+pi) & halo_elev<(pi-cuts_elev(i,1)) & halo_elev>(pi-cuts_elev(i,2));
%              cut = or(cut, cut2);
%         end
%         halo_radial=halo_radial(~cut);
%         halo_azm=halo_azm(~cut);
%         halo_elev = halo_elev(~cut);        
%      end
%     end
%     
%     %launch into the nested loop of binning
%     for m=[1:size(bin_pairs_azm,1)]
%         for p=[1:size(bin_pairs_elev,1)]
%             %count the number of atoms in bin m,p
%         
%             %separate clause for bins wrapping around 0,2pi to prevent lost counts
%             if bin_pairs_azm(m,1)>bin_pairs_azm(m,2)                
%                 %break bin into two parts (x1,2pi) and (0,x2) and count separately
%                 counts = sum( halo_azm>bin_pairs_azm(m,1) & halo_azm<2*pi ...
%                     & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2))...
%                     + sum( halo_azm>0 & halo_azm<bin_pairs_azm(m,2) ...
%                     & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2));           
%             else 
%             
%             counts=sum(halo_azm>bin_pairs_azm(m,1) & halo_azm<bin_pairs_azm(m,2)...
%                 & halo_elev>bin_pairs_elev(p,1) & halo_elev<bin_pairs_elev(p,2)); %mask_rad
%             end
%                                 
%             angle_counts_file{m,p}=[n,bin_centers_azm(m),bin_centers_elev(p),counts];        
%         end
%     end
%     angle_counts{n}=vertcat(angle_counts_file{:});
% end
% 
% %now we have the bins we can find the normalized number vairance
% %between them as defined in http://arxiv.org/pdf/1008.0845.pdf
% 
% num_bins = size(angle_counts{1},1);
% %all possible combinations of two bins is given by
% index_combs=nchoosek(1:num_bins,2);
% if isverbose
%     disp(['number of combinations of bins ',int2str(size(index_combs,1))]);
% end
% 
% %exclude combinations if bins overlap (these will have artificially low variance)
% exclude_overlapping_bins = 1;
% if exclude_overlapping_bins
%     if bin_width_azm > range(range_azm)/steps_azm || bin_width_elev > range(range_elev)/steps_elev 
%         if isverbose
%             disp(['deleting overlapping bins'])
%         end
%         index_combs = delete_overlapping_bins(index_combs, steps_azm, steps_elev, bin_width_azm, bin_width_elev, range_azm, range_elev,isverbose);
%     end
% end
% 
% %exclude combinations if any bin overlaps with removed parts of halo
% %find bins from bin_pair_azm that overlap with cuts, work out what number
% %this corresponds to and get rid of those
% %only works for azm cuts (doesn't distinguish elevation)
% exclude_cut_bins=1;
% if exclude_cut_bins
%     cut_bins_azm=[];
%     for m=1:size(bin_pairs_azm,1)
%         if bin_pairs_azm(m,1)>bin_pairs_azm(m,2)
%             interval1=fixed.Interval(bin_pairs_azm(m,1),2*pi,'()');
%             interval2=fixed.Interval(0,bin_pairs_azm(m,2),'()');
%             bin_range=union(interval1,interval2);
%         else
%         bin_range=fixed.Interval(bin_pairs_azm(m,1),bin_pairs_azm(m,2),'()'); %work out what to do with 2pi
%         end
%         %check if it overlaps with any cuts (including opposite to a cut!)
%         for i=1:size(cuts_azm,1)
%             cut_range=fixed.Interval(cuts_azm(i,1),cuts_azm(i,2),'()');
%             cut_range_opp = fixed.Interval(wrapTo2Pi(cuts_azm(i,1)+pi),wrapTo2Pi(cuts_azm(i,2)+pi),'()');
%             if any(overlaps(bin_range,cut_range)) || any(overlaps(bin_range,cut_range_opp))
%                 cut_bins_azm = [cut_bins_azm m];
% 
%                 %if it does overlap with a cut, delete bin combinations involving this bin
%                 if m==steps_azm
%                     cut_mask = any(mod(index_combs(:,:),steps_azm)==0,2) ;
%                     index_combs = index_combs(~cut_mask,:);
%                     %disp(num2str(size(index_combs,1)))
%                 else
%                     cut_mask = any(mod(index_combs(:,:),steps_azm)==m,2) ;
%                     index_combs = index_combs(~cut_mask,:);
%                     %disp(num2str(size(index_combs,1)))
%                 end
% 
%             end
%         end
% 
%     end 
% end
% 
% if isverbose
%     disp(['number of combinations of bins (cuts removed) ',num2str(size(index_combs,1))])
% end
% 
% %array to store number of counts for each bin combination for each shot
% %dimensions: number of shots, number of bin combinations, 2
% bin_pairs=zeros(size(halo_centered_cells,1),size(index_combs,1),2);
% 
% %first we rearange the cell output into an array
% %there may be a better way to do this other than a loop but i cant find it
% %is the cell array actually used for anything or can we just
% %put stuff straight into angle_counts_mat as they are calculated ?
% if isverbose
% disp('Finding Norm Var ')
% end
% 
% num_shots = size(angle_counts,2); 
% angle_counts_mat=zeros(num_bins,size(angle_counts{1},2),num_shots);
% for n=1:num_shots
%     angle_counts_mat(:,:,n)=angle_counts{n};
% end
% %angle_counts_mat [bin number, angle_counts_info, shot]
% %to get counts in bin x, shot n, angle_counts_mat(x,4,n)
% 
% %this selects the appropriate bins so we have a matix of filex*index in
% %comb list*(bin1,bin2)
% bin_pairs=[angle_counts_mat(index_combs(:,1),4,:), angle_counts_mat(index_combs(:,2),4,:)]; %[counts in bin1, counts in bin2]
% %permute these get desired format: first index = shot, second = bin comb.,
% %third index = which of the two bins in the combination
% bin_pairs=permute(bin_pairs,[3 1 2]);
% 
% %select only data from shots with number of counts
% if window_counts
%     count_mask = total_counts>count_lims(1) & total_counts<count_lims(2);
%     bin_pairs=bin_pairs(count_mask,:,:);
%     if isverbose
%         disp(['number of good shots: ',num2str(sum(count_mask))])
%     end
%     
%     if plot_count_distribution
%         stfig('Number of counts in each shot');
%         clf
%         hist(total_counts,2.5:5:max(total_counts));
%     end
% end
% 
% 
% %calculate normalised variance
% diffs=(bin_pairs(:,:,1)-bin_pairs(:,:,2));
% diffs2=(diffs).^2;
% norm_var=(mean(diffs2,1,'omitnan')-mean(diffs,1, 'omitnan').^2)./(mean(bin_pairs(:,:,1),'omitnan')+mean(bin_pairs(:,:,2), 'omitnan'));
% 
% %calculate the angle between the bin_centers 
% %https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
% %https://math.stackexchange.com/questions/243142/what-is-the-general-formula-for-calculating-dot-and-cross-products-in-spherical
% %r1r2(sin?1sin?2cos(?1??2)+cos?1cos?2)
% %acos(sin(elev1)sin(elev2)cos(azm1?azm2)+cos(elev1)cos(elev2))
% angle_pairs=zeros(1,size(index_combs,1));
% for n=[1:size(index_combs,1)]
%     azm1=angle_counts{1}(index_combs(n,1),2);
%     elev1=angle_counts{1}(index_combs(n,1),3);
%     azm2=angle_counts{1}(index_combs(n,2),2);
%     elev2=angle_counts{1}(index_combs(n,2),3);
%     
%     angle_pairs(n)= real(acos(sin(elev1)*sin(elev2)*cos(azm1-azm2)+cos(elev1)*cos(elev2)));
%     %only take real part to avoid error where acos returns a tiny imaginary component
% end
% 
% %sort through the angle and norm var values, average the norm_var values that have the same angle_pairs
% angle_tol= 0.0001;% 0.0001;
% uniq_angles=uniquetol(angle_pairs,angle_tol);
% angle_var_sd=zeros(size(uniq_angles,2),3);
% %here i also remove zeros
% for n=1:size(uniq_angles,2)
%     matching=norm_var((abs(angle_pairs-uniq_angles(n))<angle_tol) & norm_var~=0);
%     angle_var_sd(n,:)=[uniq_angles(n)/pi, mean(matching,'omitnan'),std(matching,'omitnan')/sqrt(size(matching,2))]; %omitnan
% end