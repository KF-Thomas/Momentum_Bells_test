%% AOM TRANSFER LINEARITY
% data taken around july 2020
% data = [30, 0.22, 0.056;
%     50, 0.57, 0.15;
%     100, 2.2, 0.5;
%     150, 4.85, 1.25;
%     200, 8.2, 2.2;
%     250, 12.7, 3.3;
%     300, 18.2, 4.8;
%     400, 29.8, 8.2;
%     500, 43.3, 12.6;
%     600, 57.5, 17];

% %data taken 2020-11-30
% data = [100, 0.53, 2.20;
%     200, 2.06, 8.62;
%     300, 4.63, 19.5;
%     400, 8.03, 31.3;
%     500, 12.03, 49.2;%
%     600, 16.95, 59.2
%     700, 22.2,  72];

%data taken 2021-02-01
data = [100, 0.51, 2.04
    200, 2.06, 8.08;
    300, 4.61, 17.2;
    400, 7.9, 29.1;
    500, 11.8, 41.4;%
    600, 16, 54;
    700, 21.3,  66.5;
    800, 24, 75.5;
    900, 26.8, 81.7
    1000, 28.2, 85.6
    1100, 31, 86.2
    1200, 31, 86.7];
%taken 2023-03-4
data = [100,2.56,2.1
    200,9.54,7.8
    300,20.9,16.8
    400,35.6,27.9
    500,53,42
    600,70.5,58.2
    700,84.8,72.8
    800,94.2,87.3
    900,103.5,93.7
    1000,109.7,105];

stfig('AOM Transfer Efficency');
% clf
plot(data(:,1),data(:,2),'o')
hold on
plot(data(:,1),data(:,3),'o')

fit_mask = data(:,1)<600;

modelfun = @(b,x) b(1).*x(:,1).^b(2);
modelfun_sin = @(b,x) b(1).*sin(x(:,1).*b(2)).^2;
fit_1 = fitnlm(data(fit_mask,1),data(fit_mask,2),modelfun,[1.4e-3,1.6]);
fit_2 = fitnlm(data(fit_mask,1),data(fit_mask,3),modelfun,[1.6e-4,1.8]);
fit_1_sin = fitnlm(data(:,1),data(:,2),modelfun_sin,[100,1/600]);
fit_2_sin = fitnlm(data(:,1),data(:,3),modelfun_sin,[100,1/600]);

modelfun_in = @(b,x) b(1).*x(:,1).^b(2);
fit_1in = fitnlm(data(fit_mask,2),data(fit_mask,1),modelfun_in,[1,1]);
fit_2ih = fitnlm(data(fit_mask,3),data(fit_mask,1),modelfun_in,[1,1]);

xx=linspace(100,1300,1e3)';
parms1 = fit_1.Coefficients.Estimate;
params1in = fit_1in.Coefficients.Estimate;% Get Estimated Parameters
plot(xx, modelfun(parms1,xx))
parms2 = fit_2.Coefficients.Estimate;                      % Get Estimated Parameters
plot(xx, modelfun(parms2,xx))

parms1_sin = fit_1_sin.Coefficients.Estimate;
plot(xx, modelfun_sin(parms1_sin,xx))
parms2_sin = fit_2_sin.Coefficients.Estimate;                      % Get Estimated Parameters
plot(xx, modelfun_sin(parms2_sin,xx))
legend('Bragg 1','Bragg 2')
xlabel('voltage amplitude peak to peak (mV)')
ylabel('Output power (mW)')