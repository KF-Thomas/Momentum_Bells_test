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

%data taken 2020-11-30
data = [100, 0.53, 2.20;
    200, 2.06, 8.62;
    300, 4.63, 19.5;
    400, 8.03, 31.3;
    500, 12.03, 49.2;
    600, 16.95, 59.2
    700, 22.2,  72];

stfig('AOM Transfer Efficency');
clf
plot(data(:,1),data(:,2),'o')
hold on
plot(data(:,1),data(:,3),'o')

modelfun = @(b,x) b(1).*x(:,1).^b(2);
% modelfun = @(b,x) b(1).*sin(x(:,1).*b(2)).^2;
fit_1 = fitnlm(data(:,1),data(:,2),modelfun,[1.4e-3,1.6]);
fit_2 = fitnlm(data(:,1),data(:,3),modelfun,[1.6e-4,1.8]);
% fit_1 = fitnlm(data(:,1),data(:,2),modelfun,[100,1/600]);
% fit_2 = fitnlm(data(:,1),data(:,3),modelfun,[100,1/600]);

modelfun_in = @(b,x) b(1).*x(:,1).^b(2);
fit_1in = fitnlm(data(:,2),data(:,1),modelfun_in,[1,1]);
fit_2ih = fitnlm(data(:,3),data(:,1),modelfun_in,[1,1]);

parms1 = fit_1.Coefficients.Estimate;
params1in = fit_1in.Coefficients.Estimate;% Get Estimated Parameters
plot(data(:,1), modelfun(parms1,data(:,1)))
parms2 = fit_2.Coefficients.Estimate;                      % Get Estimated Parameters
plot(data(:,1), modelfun(parms2,data(:,1)))
legend('Bragg 1','Bragg 2')
xlabel('voltage amplitude peak to peak (mV)')
ylabel('Output power (mW)')