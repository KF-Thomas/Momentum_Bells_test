%% Interfering Halos Simulation
%% Experimental Parameters

num_shots = 1000; %number of shots to take

qe = 1.0;%0.3;%0.08; %quantum efficeny

halo_rad = 0.063897; %halo radius of halo m/s

corr_len = 0.005; % corrleation length in m/s

mode_occupancy = 0.3;

phase = 0.2; %relative phase of eam splitters in radians

%% Derived parameters

lambda = sqrt(mode_occupancy/(1+mode_occupancy));

delta_ang = 2*asin(corr_len/halo_rad);

theta = 0:delta_ang:(pi-delta_ang);

phi = -pi/3:delta_ang:pi/3;

total_modes = size(theta,2)*size(phi,2);

%% Output parameters
top_halo.counts_vel = cell(num_shots,1);
btm_halo.counts_vel = cell(num_shots,1);
%%
sample_dist = sample_distribution(phase,lambda);
for shot_num = 1:num_shots
seeds = rand(total_modes);
states = zeros(total_modes,4);
ii = 1;
for itheta = 1:length(theta)
    ctheta = theta(itheta);
    for jphi = 1:length(phi)
        cphi = phi(jphi);
seed = find(seeds(ii)<sample_dist.prob_dist(2:end) & seeds(ii)>sample_dist.prob_dist(1:end-1));
states(ii,:) = sample_dist.state_vec(seed,:);
sgn = [1,-1,1,-1];
shft = [0,pi,0,pi];
zshift = halo_rad.*[1,1,-1,-1];
for p = 1:4
    N = states(ii,p);
    if N>0
        for jj=1:N
        x = halo_rad*cos(ctheta+shft(p))*sin(cphi*sgn(p))+randn().*corr_len;
        y = halo_rad*sin(ctheta+shft(p))*sin(cphi*sgn(p))+randn().*corr_len;
        z = halo_rad*cos(cphi*sgn(p))+randn().*corr_len;%+zshift(p);
        point = [z x y];
        %correlation noise
        %radius noise
        
        %qe
        if rand()<qe
            if p<3
                top_halo.counts_vel{shot_num} = [top_halo.counts_vel{shot_num};point];
            else
                btm_halo.counts_vel{shot_num} = [btm_halo.counts_vel{shot_num};point];
            end
        end
        end
    end
end
ii = ii +1;
    end
end
end

%%
%% calculate the global correlation functions around the halos
global_corrs_opts.plots = true;
global_corrs_opts.fit = true;
global_corrs_opts.calc_err = false;

corrs = global_corrs(top_halo,bottom_halo,global_corrs_opts);