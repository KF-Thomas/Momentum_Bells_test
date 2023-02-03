function squeezing_zones_filtered_plot(halo_centered_cells, szcf_out)

%{

szcf_out = squeezing_zones_calc_filtered(halo{1}.counts_vel')
squeezing_zones_filtered_plot(halo{1}.counts_vel', szcf_out);

szpf_out = squeezing_zones_plot_filtered(halo{1}.counts_vel')

Nz_test = [(2:2:50) 60:10:360]';
Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
squeezing_zones_plot(Nz_results)
%}

% Nz_test = [(2:2:50) 60:10:180]';

% if ~exist('random_throw_away_perc_list', 'var')
%     random_throw_away_perc_list = 0:0.1:0.9;
% %     random_throw_away_perc_list = [0 0.5];
% end 

if ~exist("szcf_out", "var")
    szcf_out = squeezing_zones_calc_filtered(halo_centered_cells);
end 

% rtapl = random_throw_away_perc_list;
% rtapl_counts = numel(rtapl);
rtapl_counts = size(szcf_out,2);


figure(220); clf(220);
figure(220);

colors = hsv(rtapl_counts);



for irtapl = 1:rtapl_counts
    rtap = szcf_out{irtapl}{1};
    co = colors(irtapl,:);
    
%     fprintf("Calculating " + num2str(irtapl) + " of " + num2str(rtapl_counts) ...
%         + ", throwing away " + num2str(100*(1-rtap), 2) + "% of data " + '\n\n');

    label_prefix = num2str((1-rtap)*100,3) + "$\%$ data - ";
    
    squeezing_zones_out = szcf_out{irtapl}{2};

    Nz_test = squeezing_zones_out(:,1);
    Nz_mean_corr = squeezing_zones_out(:,2);
    Nz_mean_corr_std = squeezing_zones_out(:,3);
    Nz_mean_unco = squeezing_zones_out(:,4);
    Nz_mean_unco_std = squeezing_zones_out(:,5);

    [Nz_mean_corr_mean, Nz_mean_corr_mean_ind] = min(Nz_mean_corr);
%     fprintf("Minimum V = "+num2str(Nz_mean_corr_mean) + " at Nz = " + num2str(Nz_test(Nz_mean_corr_mean_ind)));
    
    % plot(lin_model,'color','red');hold on;
    % plot(lin_model.);hold on;
    
%     errorbar(Nz_test, Nz_mean_corr, Nz_mean_corr_std, '.','CapSize',0,'DisplayName',label_prefix+'correlated'); hold on;
%     errorbar(Nz_test, Nz_mean_unco, Nz_mean_unco_std, '.','CapSize',0,'DisplayName',label_prefix+'uncorrelated'); hold on;
    errorbar(Nz_test, Nz_mean_corr, Nz_mean_corr_std,'o', ...
        'CapSize',0,'Color',[co 1.0],'HandleVisibility','off','MarkerFaceColor',co,'MarkerSize',5); hold on;
    errorbar(Nz_test, Nz_mean_unco, Nz_mean_unco_std,'.', ...
        'CapSize',0,'Color',[co 0.1],'HandleVisibility','off'); hold on;
%     set([eb.Bar, eb.Line], 'ColorType', 'truecoloralpha', 'ColorData', [eb.Line.ColorData(1:3); 255*0.2]);
%     set(eb.Marker, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [eb.Cap.EdgeColorData(1:3); 0.1]);
%     set(eb, 'MarkerFa, [co 0.1]);
%     eb.Color(4) = 0.1;

    lin_model = fitlm(Nz_test, Nz_mean_corr,"linear", "Weights",1./Nz_mean_corr_std.^2);
%     disp(lin_model);
%     [p, s] = polyfit(Nz_test, Nz_mean_corr,1);
%     [yfit, dy] = polyconf(p, Nz_test, s, 'predopt','curve');
    % [yfit, dy] = polyconf(lin_model.Coefficients.Estimate, Nz_test, lin_model.Coefficients.SE, 'predopt','curve');
    
    % line(Nz_test,yfit,'color','b','LineWidth',2);
    % line(Nz_test,yfit-dy,'color','r','linestyle',':');
    % line(Nz_test,yfit+dy,'color','r','linestyle',':');
    
    nz_lin = linspace(Nz_test(1), Nz_test(end),1000)';
    [y_mod, y_mod_ci] = predict(lin_model, nz_lin);
    plot(nz_lin, y_mod, "DisplayName", label_prefix+"Correlated Linear Fit", 'Color',[co 1.0]); hold on;
%     plot(nz_lin, y_mod_ci, '--', "Color", "r", "DisplayName", "Linear Fit CI");
    % xfill = [nz_lin, fliplr(nz_lin)];
    % fill(xfill, [y_mod_ci(:,2), fliplr(y_mod_ci(:,1))], 'red');
    % fill(nz_lin, y_mod_ci, 'red');
    fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], co, ...
        'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Linear Fit $2\sigma$ CI")

end 
hold off;
% legend('correlated', 'uncorrelated');
legend();

xlabel("number of zones $N_z$",'Interpreter','latex');
ylabel("mean of normalised variance $\langle V \rangle$",'Interpreter','latex');

shg

end





