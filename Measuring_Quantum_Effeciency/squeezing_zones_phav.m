function out=squeezing_zones_phav(halo_centered_cells, Npz_test, random_throw_away_perc)
    
if ~exist('Npz_test', 'var')
    warning("Npz_test default value is used")
    Npz_test = [(2:2:50) 60:10:180]'; % default
    % Nzp_results = squeezing_zones_phav(halo{1}.counts_vel', Nz_test);
    % Nzp_results = squeezing_zones_phav(halo{1}.counts_vel');
    % squeezing_zones_plot(Nzp_results)
    %
end
if ~exist('random_throw_away_perc', 'var')
    random_throw_away_perc = 0;
end 

Nz_test_counts = numel(Npz_test);
% Nz_results = zeros(Nz_test_counts,4);

% last_size_fprintf = 0;

Nz_results_par = cell(Nz_test_counts,1);

shift_around_list_unscaled = (0:0.05:1);
shift_tests_counts = size(shift_around_list_unscaled,2);

% Nz_results_gird = cell(Nz_test_counts, shift_tests_counts);

% parfor_progress(0);
parfor_progress(Nz_test_counts); 

parfor inz = 1:Nz_test_counts
% parfor inz = 1:Nz_test_counts
    Nz = Npz_test(inz);
    shift_around_list = shift_around_list_unscaled/Nz;

    Nz_results_gird = cell(shift_tests_counts,1);
    for isp = 1:shift_tests_counts
        shift_around = shift_around_list(isp);

        V_ij_results = squeezing_new(halo_centered_cells, false, Nz, random_throw_away_perc, shift_around);
        V_ij      = V_ij_results(:,1);
        V_ij_std  = V_ij_results(:,2);
    
        V_ij_corr = logical(V_ij_results(:,3));
        V_ij_coli = logical(V_ij_results(:,4));
    
        V_ij_both = (V_ij_corr | V_ij_coli);
    
        mean_corr     = mean(V_ij(V_ij_corr));
        mean_corr_var = mean(V_ij_std(V_ij_corr).^2);
        mean_corr_std = sqrt(mean_corr_var);
    
        mean_uncorr = mean(V_ij(~V_ij_both));
        mean_uncorr_var = mean(V_ij_std(~V_ij_both).^2);
        mean_uncorr_std = sqrt(mean_uncorr_var);
    
        mean_coli = mean(V_ij(V_ij_coli));
        mean_coli_var = mean(V_ij_std(V_ij_coli).^2);
        mean_coli_std = sqrt(mean_coli_var);
    
        %Nz_results_par{inz} = [mean_corr; mean_corr_std; mean_uncorr; mean_uncorr_std; mean_coli; mean_coli_std];
        Nz_results_gird{isp} = [mean_corr; mean_corr_std; mean_uncorr; mean_uncorr_std; mean_coli; mean_coli_std];
    end

    results_this = [Nz_results_gird{:}];
    Nz_results_par{inz} = mean(results_this,2);
    
    parfor_progress;
end 

parfor_progress(0);




Nz_results = reshape(cell2mat(Nz_results_par), 6, Nz_test_counts)';

out = [Npz_test, Nz_results];

end