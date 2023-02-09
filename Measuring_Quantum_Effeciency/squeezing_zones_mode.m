function output = squeezing_zones_mode(halo_centered_cells,plot_on, Nz_test, Mode_test, random_throw_away_perc)
    arguments
        halo_centered_cells;
        plot_on = true;
        Nz_test = [(2:2:50) 60:10:180]';
        Mode_test = [10,20,Inf];
        random_throw_away_perc = 0;
    end 
    % squeezing_zones_mode.m 
    %   scans through squeezing_zones_phav() with halo of counts using Mode_test
    % 
    % 
    % 
    %     
    %
    % 
    % squeezing_zones_mode(halo_counts_data)
    % squeezing_zones_mode(halo_counts_data, true, [(2:2:50) 60:10:180]', [10,20,Inf], 0);
    % squeezing_zones_mode(halo_counts_data, true, [(2:2:50) 60:10:180]', [5,15,Inf], 0);
    % squeezing_zones_mode(halo_counts_data, true, [(2:2:50) 60:10:180]', [5,15,25,Inf], 0);
    %
    
    mode_test_counts = numel(Mode_test);
    
    mode_test_results = cell(mode_test_counts,1);
    
    colors = hsv(mode_test_counts);
    
    mode_cut_lower = 0;
    
    if plot_on
        figure(250);clf(250);figure(250);
    end 


    for im = 1:mode_test_counts
        mode_cut_upper = Mode_test(im);

        mode_halo_ind = cellfun(@(c) size(c,1) > mode_cut_lower & size(c,1) < mode_cut_upper, halo_centered_cells);

        mode_halo = {halo_centered_cells{mode_halo_ind}};
        Nzp_results = squeezing_zones_phav(mode_halo, Nz_test, random_throw_away_perc);
        mode_test_results{im} = Nzp_results;

        %%
        if plot_on 
            squeezing_zones_out = Nzp_results;

            % Warning: below code copy/duplicate of squeezing_zones_plot
            Nz_test = squeezing_zones_out(:,1);
            Nz_mean_corr = squeezing_zones_out(:,2);
            Nz_mean_corr_std = squeezing_zones_out(:,3);
            Nz_mean_unco = squeezing_zones_out(:,4);
            Nz_mean_unco_std = squeezing_zones_out(:,5);
            Nz_mean_coli = squeezing_zones_out(:,6);
            Nz_mean_coli_std = squeezing_zones_out(:,7);

            co = colors(im,:);
            mode_str = num2str(mode_cut_lower) + " to " + num2str(mode_cut_upper);
            errorbar(Nz_test, Nz_mean_unco, Nz_mean_unco_std, 'x', ...
                'Color',co,'CapSize',0,'DisplayName',mode_str+" uncorr"); hold on;
            errorbar(Nz_test, Nz_mean_corr, Nz_mean_corr_std, '.', ...
                'Color',co,'CapSize', 0,'DisplayName',mode_str+" corr",'MarkerSize',5); hold on;

            nz_lin = linspace(1, Nz_test(end)*1.02,1000)';
            cr = 0.065*2*pi;
            F = @(arg, nn) 1 - arg(1).*(nn./cr).*((-1+exp(-cr^2./(2*(nn.*arg(2)).^2))).*sqrt(2/pi).*arg(2) + ...
               (cr./nn).*erf(cr./(sqrt(2).*nn.*arg(2))));
            mdl = fitnlm(Nz_test,Nz_mean_corr,F,[0.1, 0.065*0.05],"Weights",1./Nz_mean_corr_std.^2);
            
            [y_mod, y_mod_ci] = predict(mdl, nz_lin);
            plot(nz_lin, y_mod,"Color",co, "DisplayName", "Prob Fit");
            fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], co, ...
                'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Prob Fit $2\sigma$ CI")
            
            fprintf(mode_str + "  \t fitted η = " + num2str(mdl.Coefficients.Estimate(1)) + ...
                ' ± ' + num2str(mdl.Coefficients.SE(1)) + ...
                '\t σ = ' + num2str(mdl.Coefficients.Estimate(2)) + ...
                ' ± ' + num2str(mdl.Coefficients.SE(2)) + '\n');

        end 

        mode_cut_lower = mode_cut_upper;
    end 

    if plot_on
        xlabel("number of zones $N_z$",'Interpreter','latex');
        ylabel("mean of normalised variance $\langle V \rangle$",'Interpreter','latex');
        xlim([0,Nz_test(end)*1.02]);
        ylim([0.8, 1.2]);
        legend();
    end 


    output = mode_test_results;
   
end