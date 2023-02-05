function out=squeezing_phase(halo_centered_cells,plot_on, phase_test, zones_azm, random_throw_away_perc)

if ~exist('phase_test', 'var')
    phase_test = (0:0.25:2)'*pi; % default
%     phase_results = squeezing_phase(halo{1}.counts_vel', true, phase_test, 4,0);
    % Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
    % Nz_results = squeezing_zones(halo{1}.counts_vel',true);
end
if ~exist('random_throw_away_perc', 'var')
    random_throw_away_perc = 0;
end 
if ~exist('zones_azm', 'var')
    zones_azm = 4;
end 

phase_test_counts = numel(phase_test);
% phase_results = zeros(phase_test_counts,4);

% last_size_fprintf = 0;
% fprintf("Calculating squeezing_phase ... ");
phase_results_par = cell(phase_test_counts,1);
parfor_progress(0);
fprintf('     \n');
parfor_progress(phase_test_counts); 
parfor inz = 1:phase_test_counts 
    phase = phase_test(inz);

%     if plot_on
%         fprintf(repmat('\b', 1, last_size_fprintf));
%         last_size_fprintf = fprintf("calculating " + num2str(inz) + " out of " + num2str(phase_test_counts) + ...
%             ", phase = " + num2str(phase));
%     end 

    %     [V_ij,  V_ij_std, V_ij_corr] = squeezing_new(halo_centered_cells, false, Nz_test);
    V_ij_results = squeezing_new(halo_centered_cells, false, zones_azm, random_throw_away_perc, phase);
    V_ij      = V_ij_results(:,1);
    V_ij_std  = V_ij_results(:,2);
    V_ij_corr = logical(V_ij_results(:,3));

    mean_corr     = mean(V_ij(V_ij_corr));
    mean_corr_var = mean(V_ij_std(V_ij_corr).^2);
    mean_corr_std = sqrt(mean_corr_var);
    
    mean_uncorr     = mean(V_ij(~V_ij_corr));
    mean_uncorr_var = mean(V_ij_std(~V_ij_corr).^2);
    mean_uncorr_std = sqrt(mean_uncorr_var);

%     phase_results(inz,1) = mean_corr;
%     phase_results(inz,2) = mean_corr_std;
%     phase_results(inz,3) = mean_uncorr;
%     phase_results(inz,4) = mean_uncorr_std;
%     phase_results(inz,:) = [];
    phase_results_par{inz} = [mean_corr; mean_corr_std; mean_uncorr; mean_uncorr_std];
    parfor_progress;
end 
parfor_progress(0);

phase_results = reshape(cell2mat(phase_results_par), 4,phase_test_counts)';
% fprintf("Done"+'\n');


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

out = [phase_test, phase_results];

end

