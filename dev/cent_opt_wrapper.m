% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,1,[-0.15,0.15]);
% BestSolHistory=differential_evolution(temp_func,3,[0          0     0],[-4,-4,-4],[4,4,4])

% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,0,[-0.25,0.25]);
% BestSolHistory=differential_evolution(temp_func,3,[0          0     0],[-4,-4,-4],[4,4,4])


temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,x,0,[-0.6,0.6]);
BestSolHistory=differential_evolution(temp_func,3,[0          0     0],[-4,-4,-4],[4,4,4])

% 
% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,[x 0 0],2,[-0.15,0.15]);
% BestSolHistory=differential_evolution(temp_func,6,[0 0 0 0 0 0],[-4,-4,-4,-4,-4,-4],[4,4,4,4,4,4])

% temp_func = @(x) -halo_cent_opt(data_masked_halo,bec_masked_halo,[0.44264      3.4344      1.0304 1.8233     0.19556     0.64581 x],2);
% BestSolHistory=differential_evolution(temp_func,2,[0 0],[-0.5,-0.5],[0.5,0.5])