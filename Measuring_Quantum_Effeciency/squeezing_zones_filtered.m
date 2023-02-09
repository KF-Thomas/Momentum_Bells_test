function output=squeezing_zones_filtered(halo_centered_cells, random_throw_away_perc_list, Nz_test, shift_around)
% squeezing_zones_filtered.m
%   throw away different amount of data with amount config in random_throw_away_perc_list
% 
% 
% 
%{

szf_out = squeezing_zones_filtered(halo{1}.counts_vel')

Nz_test = [(2:2:50) 60:10:360]';
Nz_results = squeezing_zones(halo{1}.counts_vel',true, Nz_test);
squeezing_zones_plot(Nz_results)
%}

% Nz_test = [(2:2:50) 60:10:180]';

if ~exist('random_throw_away_perc_list', 'var')
    random_throw_away_perc_list = 0:0.2:0.9;
%     random_throw_away_perc_list = [0 0.5];
end 
if ~exist('Nz_test', 'var')
    Nz_test = [(2:2:50) 60:10:180]';
end 
if ~exist('shift_around','var')
    shift_around = 0;
end 

rtapl = random_throw_away_perc_list;
rtapl_counts = numel(rtapl);

output = {};

% clf(201);
% figure(201);
% rtap = 0;
% fprintf('\n');

parfor_progress(rtapl_counts); 
parfor irtapl = 1:rtapl_counts
    rtap = rtapl(irtapl);

%     fprintf("Calculating " + num2str(irtapl) + " of " + num2str(rtapl_counts) ...
%         + ", throwing away " + num2str(100*(1-rtap), 2) + "% of data " + '\n\n');
    
    output{irtapl} = {rtap, squeezing_zones(halo_centered_cells,true, Nz_test, rtap, shift_around)};
    parfor_progress;
end 
parfor_progress(0);


end