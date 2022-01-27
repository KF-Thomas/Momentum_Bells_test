y12=[]
for ii = 1:6
% y12(ii)=out_corrs{ii}.g12.norm_g2.fitted_g2peak;
% y34(ii)=out_corrs{ii}.g34.norm_g2.fitted_g2peak;
% y14(ii)=out_corrs{ii}.g14.norm_g2.fitted_g2peak;
% y23(ii)=out_corrs{ii}.g23.norm_g2.fitted_g2peak;

% y12(ii)=out_corrs{ii}.g12.in_shot_corr.rad_corr_density(1);
% y34(ii)=out_corrs{ii}.g34.in_shot_corr.rad_corr_density(1);
% y14(ii)=out_corrs{ii}.g14.in_shot_corr.rad_corr_density(1);
% y23(ii)=out_corrs{ii}.g23.in_shot_corr.rad_corr_density(1);

% y12(ii)=out_corrs{ii}.g12.norm_g2.g2_amp(2);
% y34(ii)=out_corrs{ii}.g34.norm_g2.g2_amp(2);
% y14(ii)=out_corrs{ii}.g14.norm_g2.g2_amp(2);
% y23(ii)=out_corrs{ii}.g23.norm_g2.g2_amp(2);

y12(ii)=mean(out_corrs{ii}.g12.norm_g2.g2_amp(1:3));
y34(ii)=mean(out_corrs{ii}.g34.norm_g2.g2_amp(1:3));
y14(ii)=mean(out_corrs{ii}.g14.norm_g2.g2_amp(1:3));
y23(ii)=mean(out_corrs{ii}.g23.norm_g2.g2_amp(1:3));

end

for ii = 1:6
    Ntop(ii) = mean(out_data{ii}.top_halo.num_counts);
end