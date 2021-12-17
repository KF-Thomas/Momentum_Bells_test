function [xg2,k_vec,weight] = g2_cumulative_avg(corrs,dk)
%returns a vector of g2 amplitudes for varrying integration bin size
centers = corrs.in_shot_corr.x_centers;
cen_mask = (centers == 0);
cen_indx = find(cen_mask);

raw_corr = corrs.in_shot_corr.one_d_corr_density(cen_indx:end);
raw_unc = corrs.in_shot_corr.corr_unc(cen_indx:end);
norm_corr = corrs.between_shot_corr.one_d_corr_density(cen_indx:end);
   %desired_result=cumsum(a)./(1:numel(a))
   
   
   
    numerator = cumsum(raw_corr);
    denominator = cumsum(norm_corr);
    xg2 = numerator./denominator;
    
    weight = sqrt(cumsum(raw_unc.^2))./(1:numel(raw_unc))';
    weight = weight./denominator;
    
    k_vec = cumsum(centers(cen_indx:end))./(1:numel(centers(cen_indx:end)))';
    dk_vec = (1:numel(centers(cen_indx:end))).*dk;
    
end