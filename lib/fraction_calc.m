function out = fraction_calc(data,frac_opts)
transfer_state = frac_opts.transfer_state; %indicates if we are calibrating magnetic or momentum transfer

%% run some checks
% atoms number
% laser maybe?
num_check = data.num_counts>frac_opts.num_lim;
is_shot_good = num_check;% & tag_mask;
data_masked = struct_mask(data,is_shot_good);
%% Count up in the different time bins
% time bounds we care about
if strcmp(transfer_state,'mag')
%     t_bounds = {[1.752,1.763],[1.763,1.776],[1.776,1.787],[1.75,1.79]};%time bounds for current tight trap settings
%     t_bounds = {[2.15,2.1632],[2.1632,2.175],[2.175,2.19],[2.1,2.25]};
    t_bounds = {[3.86,3.8725],[3.8725,3.881],[3.881,3.895],[3.75,4]}; %time bounds for the different magnetic states mj=+,1,0,-1 respectively
elseif strcmp(transfer_state,'momentum')
 %   t_bounds = {[3.8598,3.871],[3.871,3.8844],[3.8844,3.8972],[3.75,4]};
%     t_bounds = {[2.134,2.148],[2.148,2.161],[2.161,2.18],[2.13,2.2]};%
% t_bounds = {[1.893,1.905],[1.88,1.89],[1.865,1.877],[1.865,1.905]};%
%  t_bounds = {[1.08,1.09],[1.065,1.075],[1.04,1.05],[1.035,1.095]};%
%   t_bounds = {[2.095,2.105],[2.08,2.09],[2.065,2.075],[2.06,2.11]};% 
% t_bounds = {[2.262,2.268],[2.245,2.25],[2.253,2.259],[2.24,2.27]}; % Rabi 2023 Sep 19?
% t_bounds = {[2.30,2.31],[2.268,2.278],[2.33,2.24],[2.24,2.34]}; % Raman 2023 Sep 20?
% t_bounds = {[3.295,3.305],[3.265,3.275],[3.33,3.24],[3.24,3.34]}; % Raman 2023 Sep 20? DLD trigger super weird
t_bounds = {[2.185,2.225],[2.15,2.18],[2.3,2.32],[2.13,2.33]};
%    t_bounds = {[1.735,1.75],[1.75,1.763],[1.763,1.776],[1.73,1.779]};%
%     t_bounds = {[1.735,1.75],[1.753,1.759],[1.7665,1.772],[1.73,1.779]};%
%     t_bounds = {[3.848,3.8598],[3.8598,3.871],[3.871,3.8844],[3.815,4]};%time bounds for the different momentum states k=+,1,0,-1 respectively
end
% 
% t_bounds = {[3.8274,3.8611],[3.8611,3.9],[3.9,3.9467]}; %time bounds for the different magnetic states mj=+,1,0,-1 respectively

num_shots = length(data_masked.shot_num);
Ns = zeros(num_shots,length(t_bounds));
for shot_idx = 1:num_shots
    for t_idx = 1:length(t_bounds)
        lims = [t_bounds{t_idx};frac_opts.bounds];
        this_txy = data_masked.counts_txy{shot_idx};
        trim_txy = masktxy_square(this_txy, lims);
        Ns(shot_idx,t_idx) = size(trim_txy,1);
    end
end
Ntotals = sum(Ns(:,1:3),2);
out.fracs = Ns./Ntotals;
out.fracs(:,4) = Ns(:,4)./sum(Ns(:,1:2),2);
out.shot_num = data_masked.shot_num;
out.Ntotal = Ntotals;
out.Ns = Ns;
end