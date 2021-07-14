function bec = halo_cent_correction(bec,opts_cent_correction)
%% Centering correction from expected halo radius (currently only looks at time)

hebec_constants
vr = -0.130159/2;
g = -const.g0;
d = const.fall_distance;
tf = sqrt(-2*d/g);
delt = 1.769-tf;%3.8772-tf;%2.1749-tf;%

opposite_pol = @(t,s) (2*d+t.^2.*g+4.*t.*vr*s-sqrt(-8*t.^2.*d.*g+(2*d+t.^2.*g+4*t.*vr*s).^2))./(2*t.*g);

top_t_predicted_1 = opposite_pol(bec.centre_mid(:,1)-delt,1);
top_t_predicted_2 = opposite_pol(bec.centre_btm(:,1)-delt,2);
top_t = bec.centre_top(:,1)-delt;
top_cen_predictions = [top_t_predicted_1, top_t_predicted_2, top_t];

mid_t_predicted_1 = opposite_pol(bec.centre_top(:,1)-delt,-1);
mid_t_predicted_2 = opposite_pol(bec.centre_btm(:,1)-delt,1);
mid_t = bec.centre_mid(:,1)-delt;
mid_cen_predictions = [mid_t_predicted_1, mid_t_predicted_2, mid_t];

btm_t_predicted_1 = opposite_pol(bec.centre_mid(:,1)-delt,-1);
btm_t_predicted_2 = opposite_pol(bec.centre_top(:,1)-delt,-2);
btm_t = bec.centre_btm(:,1)-delt;
btm_cen_predictions = [btm_t_predicted_1, btm_t_predicted_2, btm_t];

top_outlier_indx = isoutlier(top_cen_predictions);
mid_outlier_indx = isoutlier(mid_cen_predictions);
btm_outlier_indx = isoutlier(btm_cen_predictions);

for jj = 1:length(bec.centre_mid(:,1))
    
    bec.centre_top(jj,1) = nanmean(top_cen_predictions(jj,~top_outlier_indx(jj,:)))+delt;
    bec.centre_mid(jj,1) = nanmean(mid_cen_predictions(jj,~mid_outlier_indx(jj,:)))+delt;
    bec.centre_btm(jj,1) = nanmean(btm_cen_predictions(jj,~btm_outlier_indx(jj,:)))+delt;
    
    xy_top_predict = -(bec.centre_top(jj,1)-delt).*(bec.centre_btm(jj,2:3)./(bec.centre_btm(jj,1)-delt)...
        -2.*bec.centre_mid(jj,2:3)./(bec.centre_mid(jj,1)-delt));
    
    nan_indx = isnan(bec.centre_top(jj,2:3));
    if sum(nan_indx)>0
        bec.centre_top(jj,logical([0 nan_indx])) = xy_top_predict(nan_indx);
    end
    
    xy_mid_predict = (bec.centre_mid(jj,1)-delt).*(bec.centre_btm(jj,2:3)./(bec.centre_btm(jj,1)-delt)...
        +bec.centre_top(jj,2:3)./(bec.centre_top(jj,1)-delt));
    nan_indx = isnan(bec.centre_mid(jj,2:3));
    if sum(nan_indx)>0
        bec.centre_mid(jj,logical([0 nan_indx])) = xy_mid_predict(nan_indx);
    end
    
    xy_btm_predict = -(bec.centre_btm(jj,1)-delt).*(bec.centre_top(jj,2:3)./(bec.centre_top(jj,1)-delt)...
        -2.*bec.centre_mid(jj,2:3)./(bec.centre_mid(jj,1)-delt));
    nan_indx = isnan(bec.centre_btm(jj,2:3));
    if sum(nan_indx)>0
        bec.centre_btm(jj,logical([0 nan_indx])) = xy_btm_predict(nan_indx);
    end
    
    if all(~isnan(bec.centre_mid(jj,:)))
        bec.centre_OK_mid(jj) = 1;
    end
    if all(~isnan(bec.centre_top(jj,:)))
        bec.centre_OK_top(jj) = 1;
    end
    if all(~isnan(bec.centre_btm(jj,:)))
        bec.centre_OK_btm(jj) = 1;
    end
end
if opts_cent_correction.plots
    stfig('Bottom centering comparisson')
    clf
    plot(bec.centre_btm(:,1)-delt,'o')
    hold on
    plot(btm_t_predicted_1,'x');
    plot(btm_t_predicted_2,'.');
    plot(indxs(outlier_indx),btm_cen_predictions(outlier_indx),'rs')
    legend('btm','mid','top')
    
    stfig('Top centering comparisson')
    clf
    plot(bec.centre_top(:,1)-delt,'o')
    hold on
    plot(top_t_predicted_1,'x');
    plot(top_t_predicted_2,'s');
    
    stfig('Mid centering comparisson')
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
end
end