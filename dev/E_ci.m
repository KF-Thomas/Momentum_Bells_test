function ci = E_ci(ports,opts_E)
%calculates the confidence interval of tthe quantum correlator E for our Bell test
%input data as structure of the 4 different ports

est_fun = @(x) corr_wrap(ports,opts_E,x);
m=(1:length(ports.bottom_right.norm))';
if isfield(opts_E,'confidence_interval')
    ci = bootci(opts_E.num_samp_rep,est_fun,m,'Alpha',1-opts_E.confidence_interval);
else
    ci = bootci(opts_E.num_samp_rep,{est_fun,m},'Alpha',1-0.6826895,'Type','norm');
end

end

function out = corr_wrap(ports,opts_E,x)
ports_masked.top_left=struct_mask(ports.top_left,x);
ports_masked.top_right=struct_mask(ports.top_right,x);
ports_masked.bottom_left=struct_mask(ports.bottom_left,x);
ports_masked.bottom_right=struct_mask(ports.bottom_right,x);

opts_E.plots = 0;
opts_E.verbose = 0;
opts_E.calc_err = 0;
[E_val, corrs] = E(ports_masked,opts_E);

out = E_val;

end