function squeezing_zones_plot(squeezing_zones_out)

%{
Nz_test = [(2:2:50) 60:10:360]';
Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
squeezing_zones_plot(Nz_results)
%}

Nz_test = squeezing_zones_out(:,1);
Nz_mean_corr = squeezing_zones_out(:,2);
Nz_mean_corr_std = squeezing_zones_out(:,3);
Nz_mean_unco = squeezing_zones_out(:,4);
Nz_mean_unco_std = squeezing_zones_out(:,5);

[Nz_mean_corr_mean, Nz_mean_corr_mean_ind] = min(Nz_mean_corr);
fprintf("Minimum V = "+num2str(Nz_mean_corr_mean) + " at Nz = " + num2str(Nz_test(Nz_mean_corr_mean_ind)));

figure(210);clf(210);
figure(210);

% plot(lin_model,'color','red');hold on;
% plot(lin_model.);hold on;

errorbar(Nz_test, Nz_mean_corr, Nz_mean_corr_std, '.', 'Color','red','CapSize', 0,'DisplayName','correlated'); hold on;
errorbar(Nz_test, Nz_mean_unco, Nz_mean_unco_std, '.', 'Color','blue','CapSize',0,'DisplayName','uncorrelated'); hold on;
% legend('correlated', 'uncorrelated');
legend('AutoUpdate','off');


xlabel("number of zones $N_z$",'Interpreter','latex');
ylabel("mean of normalised variance $\langle V \rangle$",'Interpreter','latex');


lin_model = fitlm(Nz_test, Nz_mean_corr,"linear", "Weights",1./Nz_mean_corr_std.^2);
disp(lin_model);
[p, s] = polyfit(Nz_test, Nz_mean_corr,1);
[yfit, dy] = polyconf(p, Nz_test, s, 'predopt','curve');
% [yfit, dy] = polyconf(lin_model.Coefficients.Estimate, Nz_test, lin_model.Coefficients.SE, 'predopt','curve');

% line(Nz_test,yfit,'color','b','LineWidth',2);
% line(Nz_test,yfit-dy,'color','r','linestyle',':');
% line(Nz_test,yfit+dy,'color','r','linestyle',':');

nz_lin = linspace(Nz_test(1), Nz_test(end),1000)';
[y_mod, y_mod_ci] = predict(lin_model, nz_lin);
plot(nz_lin, y_mod,"Color","r", "DisplayName", "Linear Fit");
% plot(nz_lin, y_mod_ci, '--', "Color", "r", "DisplayName", "Linear Fit CI");
% xfill = [nz_lin, fliplr(nz_lin)];
% fill(xfill, [y_mod_ci(:,2), fliplr(y_mod_ci(:,1))], 'red');
% fill(nz_lin, y_mod_ci, 'red');
fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], 'red', ...
    'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Linear Fit $2\sigma$ CI")

shg

end