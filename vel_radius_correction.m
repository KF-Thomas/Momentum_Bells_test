hebec_constants
vr = -0.130159/2;
tf = 0.417;
g = -const.g0;
d = const.fall_distance;
delt = 3.8772-0.417;

opposite_pol = @(t,s) (2*d+t.^2.*g+4.*t.*vr*s-sqrt(-8*t.^2.*d.*g+(2*d+t.^2.*g+4*t.*vr*s).^2))./(2*t.*g);

top_t_predicted_1 = opposite_pol(bec.centre_mid(:,1)-delt,1);
top_t_predicted_2 = opposite_pol(bec.centre_btm(:,1)-delt,2);
stfig('Top centering comparisson')
clf
% subplot(2,1,1)
plot(bec.centre_top(:,1)-delt,'o')
hold on
plot(top_t_predicted_1,'x');
plot(top_t_predicted_2,'s');
% subplot(2,1,2)

btm_t_predicted_1 = opposite_pol(bec.centre_mid(:,1)-delt,-1);
btm_t_predicted_2 = opposite_pol(bec.centre_top(:,1)-delt,-2);
btm_cen_predictions = [btm_t_predicted_1, btm_t_predicted_2, (bec.centre_btm(:,1)-delt)];
stfig('Bottom centering comparisson')
clf
plot(bec.centre_btm(:,1)-delt,'o')
hold on
plot(btm_t_predicted_1,'x');
plot(btm_t_predicted_2,'.');
indxs = [1:length(bec.centre_mid(:,1))',1:length(bec.centre_mid(:,1))',1:length(bec.centre_mid(:,1))'];
outlier_indx = isoutlier(btm_cen_predictions);
plot(indxs(outlier_indx),btm_cen_predictions(outlier_indx),'rs')
nanmean(btm_cen_predictions(~outlier_indx(5,:)));
legend('btm','mid','top')

stfig('Mid comparisson')
mid_t_predicted_1 = opposite_pol(bec.centre_top(:,1)-delt,-1);
mid_t_predicted_2 = opposite_pol(bec.centre_btm(:,1)-delt,1);
clf
plot(mid_t_predicted_1,'o')
hold on
plot(mid_t_predicted_2,'x');
plot(bec.centre_mid(:,1)-delt,'s');
legend('top','btm','mid')

stfig('centering dif')
clf
plot(bec.centre_top(:,1)-delt-top_t_predicted_1,'o')
hold on
plot(-bec.centre_btm(:,1)+delt+btm_t_predicted_1,'x');