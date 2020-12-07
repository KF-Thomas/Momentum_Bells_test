temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,0);
BestSolHistory=differential_evolution(temp_func,3,[2.4268           0     0.23856],[-4,-4,-4],[4,4,4])

% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,2);
% BestSolHistory=differential_evolution(temp_func,6,[0 0 0 0 0 0],[-4,-4,-4,-4,-4,-4],[4,4,4,4,4,4])

% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,2);
% BestSolHistory=differential_evolution(temp_func,2,[0 0],[-0.5,-0.5],[0.5,0.5])