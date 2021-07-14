function [out, boot_out] = corr_err_raw(corr_func,corr_opts,counts)
warning('off','all')
est_fun = @(c,a,opts) corr_wrap(corr_func,c,a,opts);
if ~isfield(corr_opts,'samp_frac_lims')
    corr_opts.samp_frac_lims = [0.25,0.5];
end
if ~isfield(corr_opts,'num_samp_frac')
    corr_opts.num_samp_frac = 5;
end
if ~isfield(corr_opts,'num_samp_rep')
    corr_opts.num_samp_rep = 3e1;
end

boot_out=bootstrap_se(est_fun,counts,...
    'samp_frac_lims',corr_opts.samp_frac_lims,...
    'num_samp_frac',corr_opts.num_samp_frac,...
    'num_samp_rep',corr_opts.num_samp_rep,...
    'replace',false,...
    'use_frac_size',true,...
    'opp_arguments',{corr_opts},...
    'verbose',0);

grid_size = size(corr_opts.redges,2) - 1;
out.rawcorr.val = nansum(boot_out.sampling.mean.*boot_out.sampling.sample_size)./nansum(boot_out.sampling.sample_size);


norm_sort_dir=corr_opts.sorted_dir;
if strcmp(corr_opts.sampling_method,'basic')
    if size(counts,1)>1
        %set the number of chunks to be at least as many as the heighest
        %count number
        chunk_num = max([cellfun(@(x)size(x,1),counts(1,:)),cellfun(@(x)size(x,1),counts(2,:))]);
        counts_chunked(1,:)=chunk_data(counts(1,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
        counts_chunked(2,:)=chunk_data(counts(2,:),corr_opts.norm_samp_factor,norm_sort_dir,chunk_num);
    else
        counts_chunked=chunk_data(counts,corr_opts.norm_samp_factor,norm_sort_dir);
    end
    corr_opts.normalisation_factor = (size(counts,2)/nanmean(cellfun(@(x)size(x,1),counts)))^2;
elseif strcmp(corr_opts.sampling_method,'complete')
    counts_chunked=chunk_data_complete(counts,corr_opts.sample_proportion,norm_sort_dir);
end
normcorr=corr_func(corr_opts,counts_chunked);

if isfield(normcorr,'one_d_corr_density')
    density='one_d_corr_density';
    centers='x_centers';
elseif isfield(normcorr,'rad_corr_density')
    density='rad_corr_density';
    centers='rad_centers';
end

out.normcorr.val = normcorr.(density)';

out.corr_density.val = out.rawcorr.val./out.normcorr.val;

out.rawcorr.unc = boot_out.results.se_fun_whole_unweighted(1,1:grid_size);
out.corr_density.unc = out.rawcorr.unc./out.normcorr.val;

warning('on','all')
end

function out = corr_wrap(corr_func,counts,atten,corr_opts)
corr_opts.attenuate_counts = atten;
corr_opts.plots = false;
corr_opts.print_update = false;
rawcorr = corr_func(corr_opts,counts);

if isfield(rawcorr,'one_d_corr_density')
    density='one_d_corr_density';
    centers='x_centers';
elseif isfield(rawcorr,'rad_corr_density')
    density='rad_corr_density';
    centers='rad_centers';
end


out = [rawcorr.(density)];

end