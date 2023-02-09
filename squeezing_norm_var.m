function out=squeezing_norm_var(bin_pairs,total_counts,angle_pairs,window_counts,count_lims,plot_count_distribution,isverbose)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%calulates normalized variance for squeezing
%inputs 
    %bin_pairs ; array containing counts for each bin combination for each shot
    %total_counts ; array containing total counts in each shot
    %angle_pairs ; array containing angle between each bin pair
    %(these are outputs of squeezing_bins.m)
    
    %window_counts ; only take shots with acceptable number of total counts
    %count_lims ; acceptable range of total counts
    %plot_count_distribution ; plot histogram of total counts for each shot
    %isverbose = 1 displays some info, =0 does not. 

    %outputs
    %a cell that contains three matrices:
    %the first contains what is needed for the sqz angle plot (angle,variance,variance_se)
    %the second contains [mean_opst_bin,unc_opst_bin,mean_other_bin,unc_other_bin]
    %the third contains what is needed for the sqz_bins plot (indexed norm_var for
    %correlated/uncorrelated bins)
%--------------------------------------------------------------------------
if window_counts
    count_mask = total_counts>count_lims(1) & total_counts<count_lims(2);
    bin_pairs=bin_pairs(count_mask,:,:);
    if isverbose
        disp(['number of good shots: ',num2str(sum(count_mask))])
    end
    
    
end
if plot_count_distribution
        stfig('Number of counts in each shot');
        clf
        hist(total_counts,2.5:5:max(total_counts));
end

%calculate norm_var
diffs=(bin_pairs(:,:,1)-bin_pairs(:,:,2));
diffs2=(diffs).^2;
norm_var=(mean(diffs2,1,'omitnan')-mean(diffs,1, 'omitnan').^2)./(mean(bin_pairs(:,:,1),'omitnan')+mean(bin_pairs(:,:,2), 'omitnan'));

%sort through the angle and norm var values, average the
%norm_var values that have the same angle_pairs
angle_tol= 0.0001;% 0.0001;
uniq_angles=uniquetol(angle_pairs,angle_tol);
angle_var_sd=zeros(size(uniq_angles,2),3);
%here i also remove zeros
for n=1:size(uniq_angles,2)
    matching=norm_var((abs(angle_pairs-uniq_angles(n))<angle_tol) & norm_var~=0);
    angle_var_sd(n,:)=[uniq_angles(n)/pi, mean(matching,'omitnan'),std(matching,'omitnan')/sqrt(size(matching,2))]; %omitnan
end

%sort into opst and rest (David's old code)
norm_var_opst=norm_var((abs(angle_pairs-pi)<angle_tol) & norm_var~=0);
norm_var_rest=norm_var(~(abs(angle_pairs-pi)<angle_tol) & norm_var~=0);

%sort into opst and rest but preserve bin pair number
norm_var_pair_num = [(1:length(norm_var))' norm_var'];
mask_opst = (abs(angle_pairs-pi)<angle_tol) & norm_var~=0;

norm_var_pair_num_opst = norm_var_pair_num(mask_opst,:);
norm_var_pair_num_rest = norm_var_pair_num(~mask_opst,:);

norm_var_avg=angle_var_sd(:,1);

mean_opst_bin=angle_var_sd(abs(angle_var_sd(:,1)-1)<angle_tol,2);
unc_opst_bin=angle_var_sd(abs(angle_var_sd(:,1)-1)<angle_tol,3);
non_opst_bins=angle_var_sd(angle_var_sd(:,1)~=1,2);
mean_other_bin=mean(non_opst_bins,'omitnan'); %omitnan
min_other_bin=min(non_opst_bins);
unc_other_bin=std(non_opst_bins,'omitnan')/sqrt(length(mean_other_bin)); %omitnan
min_opst_bin=min(norm_var(angle_pairs==pi));

if  isverbose
    disp(['mean opst bin ' ,num2str(mean_opst_bin),'±',num2str(unc_opst_bin)])
    disp(['mean other bins ',num2str(mean_other_bin),'±',num2str(unc_other_bin)])
    disp(['min opst bin ',num2str(min_opst_bin)])
end

out={[angle_var_sd(:,1),angle_var_sd(:,2),angle_var_sd(:,3)],...
    [[mean_opst_bin,unc_opst_bin,min_opst_bin];[mean_other_bin,unc_other_bin,min_other_bin]],...
    {norm_var_pair_num_opst,norm_var_pair_num_rest}};
end