function out=squeezing_zones(halo_centered_cells,plot_on, Nz_test, random_throw_away_perc)

if ~exist('Nz_test', 'var')
    Nz_test = [(2:2:50) 60:10:360]'; % default
    % Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
    % Nz_results = squeezing_zones(halo{1}.counts_vel',true);
end
if ~exist('random_throw_away_perc', 'var')
    random_throw_away_perc = 0;
end 

Nz_test_counts = numel(Nz_test);
Nz_results = zeros(Nz_test_counts,4);

last_size_fprintf = 0;
for inz = 1:Nz_test_counts
    Nz = Nz_test(inz);
    if plot_on
        fprintf(repmat('\b', 1, last_size_fprintf));
        last_size_fprintf = fprintf("calculating " + num2str(inz) + " out of " + num2str(Nz_test_counts) + ...
            ", zones_azm = " + num2str(Nz));
    end 

    %     [V_ij,  V_ij_std, V_ij_corr] = squeezing_new(halo_centered_cells, false, Nz_test);
    V_ij_results = squeezing_new(halo_centered_cells, false, Nz, random_throw_away_perc);
    V_ij      = V_ij_results(:,1);
    V_ij_std  = V_ij_results(:,2);
    V_ij_corr = logical(V_ij_results(:,3));

    mean_corr     = mean(V_ij(V_ij_corr));
    mean_corr_var = mean(V_ij_std(V_ij_corr).^2);
    mean_corr_std = sqrt(mean_corr_var);
    
    mean_uncorr     = mean(V_ij(~V_ij_corr));
    mean_uncorr_var = mean(V_ij_std(~V_ij_corr).^2);
    mean_uncorr_std = sqrt(mean_uncorr_var);

    Nz_results(inz,1) = mean_corr;
    Nz_results(inz,2) = mean_corr_std;
    Nz_results(inz,3) = mean_uncorr;
    Nz_results(inz,4) = mean_uncorr_std;
end 
if plot_on 
    fprintf(repmat(' ',1, 100));
    fprintf("Done \n");
end

% if plot_on
%     figure(200);
%     errorbar(Nz_test, Nz_results(:,1), Nz_results(:,2), '.', 'Color','red','CapSize', 0); hold on;
%     errorbar(Nz_test, Nz_results(:,3), Nz_results(:,4), '.', 'Color','blue','CapSize',0); hold on;
%     disp("DEPRECATED! use squeezing_zones_plot.m")
% end 

out = [Nz_test, Nz_results];

end

