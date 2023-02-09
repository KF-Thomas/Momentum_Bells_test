function squeezing_phase_plot(squeezing_zones_out)
% plotting for results from squeezing_phase
% 
% 
% 
%{
phase_test = (0:0.01:2)'*pi; % default
phase_results = squeezing_phase(halo{1}.counts_vel', true, phase_test,4,0);
squeezing_phase_plot(phase_results)

squeezing_phase_plot(squeezing_phase(halo{1}.counts_vel', true, (0:0.01:2)'*pi,4,0))
squeezing_phase_plot(squeezing_phase(halo_counts_data, true, (0:0.01:2)'*pi,4,0))
%}

phase_test = squeezing_zones_out(:,1);
phase_mean_corr = squeezing_zones_out(:,2);
phase_mean_corr_std = squeezing_zones_out(:,3);
phase_mean_unco = squeezing_zones_out(:,4);
phase_mean_unco_std = squeezing_zones_out(:,5);

[phase_mean_corr_mean, phase_mean_corr_mean_ind] = min(phase_mean_corr);
fprintf("Minimum V = "+num2str(phase_mean_corr_mean) + " at phase = " + num2str(phase_test(phase_mean_corr_mean_ind)));

figure(211);clf(211);
figure(211);

% plot(lin_model,'color','red');hold on;
% plot(lin_model.);hold on;

errorbar(phase_test/pi, phase_mean_corr, phase_mean_corr_std, '.', 'Color','red','CapSize', 0,'DisplayName','correlated'); hold on;
errorbar(phase_test/pi, phase_mean_unco, phase_mean_unco_std, '.', 'Color','blue','CapSize',0,'DisplayName','uncorrelated'); hold on;
% legend('correlated', 'uncorrelated');
legend('AutoUpdate','off');


xlabel("Phase /$\pi$",'Interpreter','latex');
ylabel("mean of normalised variance $\langle V \rangle$",'Interpreter','latex');


lin_model = fitlm(phase_test, phase_mean_corr,"linear", "Weights",1./phase_mean_corr_std.^2);
disp(lin_model);
[p, s] = polyfit(phase_test, phase_mean_corr,1);
[yfit, dy] = polyconf(p, phase_test, s, 'predopt','curve');
% [yfit, dy] = polyconf(lin_model.Coefficients.Estimate, Nz_test, lin_model.Coefficients.SE, 'predopt','curve');

cr = [1 0.2 0.2];
mean_corr = mean(phase_mean_corr);
mean_corr_std = mean(phase_mean_corr_std);
yline(mean_corr,   '--',num2str(mean_corr,  4)+"$\pm$"+num2str(mean_corr_std,  4), ...
        'LineWidth',0.5,'Color',[cr,.9], 'Interpreter','latex'); hold on; 
rectangle('Position',[0,mean_corr-mean_corr_std   2 2*mean_corr_std], ...
        'LineStyle','none','FaceColor',[cr,.05]); hold on;

% line(Nz_test,yfit,'color','b','LineWidth',2);
% line(Nz_test,yfit-dy,'color','r','linestyle',':');
% line(Nz_test,yfit+dy,'color','r','linestyle',':');

nz_lin = linspace(phase_test(1), phase_test(end),1000)';
[y_mod, y_mod_ci] = predict(lin_model, nz_lin);
plot(nz_lin/pi, y_mod,"Color","r", "DisplayName", "Linear Fit");
% plot(nz_lin, y_mod_ci, '--', "Color", "r", "DisplayName", "Linear Fit CI");
% xfill = [nz_lin, fliplr(nz_lin)];
% fill(xfill, [y_mod_ci(:,2), fliplr(y_mod_ci(:,1))], 'red');
% fill(nz_lin, y_mod_ci, 'red');
fill([nz_lin/pi; flipud(nz_lin/pi)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], 'red', ...
    'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Linear Fit $2\sigma$ CI")

shg

end