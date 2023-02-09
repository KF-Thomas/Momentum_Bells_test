function output=sim_halo_exp_samples(radius, shots, hits, detector_noise, cut_pi_ang, cut_pi_per)
    arguments
        radius = [0.065 0.065*0.05];
        shots = 10000;
        hits = [5 0.0001];
        detector_noise = [0, 0];
        cut_pi_ang = 0;
        cut_pi_per = 0.5;
    end 

    % halo_sim = sim_halo_exp_samples();
    % squeezing_new(halo_sim.cartesian',true,10,0,0);
    % squeezing_new(halo_sim.cartesian',true,10,0.5,0);

    spherical = cell(shots,1);
    cartesian = cell(shots,1);

    for is = 1:shots

        this_shot_hits = round(normrnd(hits(1), hits(2)));
        rand_rad  = normrnd(radius(1), radius(2), this_shot_hits, 1);
%         rand_elev = 2*asin(sqrt(rand(this_shot_hits,1)));
        rand_elev = 1.8*asin(sqrt(rand(this_shot_hits,1)));
        rand_azm  = 2*pi*rand(this_shot_hits,1);
        
        rand_rad  = [rand_rad ; rand_rad];
        rand_elev = [rand_elev; pi-rand_elev];
        rand_azm  = [rand_azm ; wrapTo2Pi(rand_azm+pi)];
        
%         spherical{is} = num2cell([rand_rad, rand_azm, rand_elev]);
        spherical{is} = [rand_rad, rand_azm, rand_elev]; 

        filter_index = (rand_azm < cut_pi_ang) .* (rand(2*this_shot_hits,1) < cut_pi_per); 
%         filter_index = ((rand_azm < cut_pi_ang & rand_azm < pi) | (rand_azm < (pi+cut_pi_ang)) & rand_azm > pi) .* (rand(2*this_shot_hits,1) < cut_pi_per); 

        rand_rad  = rand_rad( ~filter_index);
        rand_elev = rand_elev(~filter_index);
        rand_azm  = rand_azm( ~filter_index);
        final_counts = size(rand_rad,1);

        x = rand_rad .* sin(rand_azm).*sin(rand_elev) + normrnd(detector_noise(1), detector_noise(2), final_counts,1);
        y = rand_rad .* cos(rand_azm).*sin(rand_elev) + normrnd(detector_noise(1), detector_noise(2), final_counts,1);
        z = rand_rad .* cos(rand_elev) + normrnd(detector_noise(1), detector_noise(2), final_counts,1);

%         cartesian{is} = num2cell([x, y, z]);
%         cartesian{is} = [x, y, z];
        
        cartesian{is} = [z, x, y];
    end

    
    
    output = struct('spherical', {spherical}, 'cartesian', {cartesian});

end

