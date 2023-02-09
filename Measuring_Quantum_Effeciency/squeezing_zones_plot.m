function squeezing_zones_plot(squeezing_zones_out, fig_id)
% squeezing_zones_plot.m 
%   Plotting results from squeezing_zones()
%
%
%
%
%
%{
Nz_test = [(2:2:50) 60:10:360]';
Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
squeezing_zones_plot(Nz_results)
%}

if ~exist("fig_id", "var")
    fig_id = 210;
end 

Nz_test = squeezing_zones_out(:,1);
Nz_mean_corr = squeezing_zones_out(:,2);
Nz_mean_corr_std = squeezing_zones_out(:,3);
Nz_mean_unco = squeezing_zones_out(:,4);
Nz_mean_unco_std = squeezing_zones_out(:,5);
Nz_mean_coli = squeezing_zones_out(:,6);
Nz_mean_coli_std = squeezing_zones_out(:,7);

[Nz_mean_corr_mean, Nz_mean_corr_mean_ind] = min(Nz_mean_corr);
fprintf("Minimum V = "+num2str(Nz_mean_corr_mean) + " at Nz = " + num2str(Nz_test(Nz_mean_corr_mean_ind)));

figure(fig_id);clf(fig_id);
figure(fig_id);

% plot(lin_model,'color','red');hold on;
% plot(lin_model.);hold on;

errorbar(Nz_test, Nz_mean_unco, Nz_mean_unco_std, '.', 'Color','blue','CapSize',0,'DisplayName','uncorrelated'); hold on;
errorbar(Nz_test, Nz_mean_coli, Nz_mean_coli_std, '.', 'Color','green','CapSize',0,'DisplayName','colinear'); hold on;
errorbar(Nz_test, Nz_mean_corr, Nz_mean_corr_std, '.', 'Color','red','CapSize', 0,'DisplayName','back to back'); hold on;
% legend('correlated', 'uncorrelated');
% legend('AutoUpdate','off');


xlabel("number of zones $N_z$",'Interpreter','latex');
ylabel("mean of normalised variance $\langle V \rangle$",'Interpreter','latex');


% lin_model = fitlm(Nz_test, Nz_mean_corr,"linear", "Weights",1./Nz_mean_corr_std.^2);
% disp(lin_model);
% [p, s] = polyfit(Nz_test, Nz_mean_corr,1);
% [yfit, dy] = polyconf(p, Nz_test, s, 'predopt','curve');
% [yfit, dy] = polyconf(lin_model.Coefficients.Estimate, Nz_test, lin_model.Coefficients.SE, 'predopt','curve');

% line(Nz_test,yfit,'color','b','LineWidth',2);
% line(Nz_test,yfit-dy,'color','r','linestyle',':');
% line(Nz_test,yfit+dy,'color','r','linestyle',':');

% nz_lin = linspace(Nz_test(1), Nz_test(end),1000)';
nz_lin = linspace(1, Nz_test(end)*1.02,1000)';
% [y_mod, y_mod_ci] = predict(lin_model, nz_lin);
% plot(nz_lin, y_mod,"Color","r", "DisplayName", "Linear Fit");
% fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], 'red', ...
%     'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Linear Fit $2\sigma$ CI")

cr = 0.065*2*pi;
% F = @(arg,xd) 1-arg(1).*(1./cr).*(xd.*(-1+exp(-cr.^2./(2*xd.^2.*arg(2)^2))).*sqrt(2./pi).*arg(2)+ ...
%     (cr./xd).*erf(cr./(sqrt(2).*xd.*arg(2))));
F = @(arg, nn) 1 - arg(1).*(nn./cr).*((-1+exp(-cr^2./(2*(nn.*arg(2)).^2))).*sqrt(2/pi).*arg(2) + ...
               (cr./nn).*erf(cr./(sqrt(2).*nn.*arg(2))));

% mdl = fitnlm(Nz_test,Nz_mean_corr,F,[0.1, 0.065*0.05],"Weights",1./Nz_mean_corr_std.^2);
mdl = fitnlm(Nz_test,Nz_mean_corr,F,[0.1, 0.065*0.05],"Weights",1./Nz_mean_corr_std.^2);

disp(mdl);

[y_mod, y_mod_ci] = predict(mdl, nz_lin);
plot(nz_lin, y_mod,"Color","r", "DisplayName", "Prob Fit");
fill([nz_lin; flipud(nz_lin)], [y_mod_ci(:,1); flipud(y_mod_ci(:,2))], 'red', ...
    'FaceAlpha',0.1, 'LineStyle','none','DisplayName',"Prob Fit $2\sigma$ CI")

% plot(nz_lin, F([0.1, 0.065*0.05], nz_lin),'-');

fprintf("At Nz = 1, fit = " + num2str(predict(mdl, 1)) + '\n');
fprintf("At Nz = 4, fit = " + num2str(predict(mdl, 4)) + '\n');
fprintf("Fitted eta = " + num2str(mdl.Coefficients.Estimate(1)) + '\n');

legend();

ylim([0.8, 1.2]);

xlim([0, Nz_test(end)*1.02]);

shg

end

















