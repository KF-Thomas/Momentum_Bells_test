%%function to sample from the interfernce of two squeeze mode states
function out =  sample_distribution(phi,lmd)
    Nlim = 5; %cu off numb
    state_vec = zeros((Nlim+1)*(Nlim+2)*(2*Nlim+3)/6,4);
    prob_dist = zeros((Nlim+1)*(Nlim+2)*(2*Nlim+3)/6+1,1);
    val = 0;
    ii = 1;
    for T = 0:Nlim
       for Y = 0:T
           for X = 0:T
               W = T-Y;
               Z = T-X;
               prob = abs(interfence_coefficent(Y,X,W,Z,phi,lmd))^2;
               if prob> eps
               val = val + prob;
               prob_dist(ii+1,:) = val;
               state_vec(ii,:) = [Y,X,W,Z];
               ii = ii + 1;
               end
           end
       end
    end
    state_vec = state_vec(1:(ii-1),:);
    prob_dist = prob_dist(1:ii,:);
    prob_dist(end,:) = 1;
    out.state_vec = state_vec;
    out.prob_dist = prob_dist;
end
