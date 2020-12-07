function [left_port, right_port] = separate_ports(data,opts_ports)
%takes in the halo in velocity space and returns the two separate ports
%indexes [z , x , y] = [3 , 1 , 2]
% split along the x-axis
left_port.unnorm = cell(size(data.counts_vel));
right_port.unnorm = cell(size(data.counts_vel));

left_port.norm = cell(size(data.counts_vel));
right_port.norm = cell(size(data.counts_vel));
for ii = 1:size(data.counts_vel,1)
    current_shot = data.counts_vel{ii};
    current_shot_norm = data.counts_vel_norm{ii};
    mask = current_shot(:,2)<0;
    mask_norm = current_shot(:,2)<0;
    left_port.unnorm{ii} = current_shot(mask,:);
    right_port.unnorm{ii} = current_shot(~mask,:);
    
    left_port.norm{ii} = current_shot_norm(mask_norm,:);
    right_port.norm{ii} = current_shot_norm(~mask_norm,:);
end

end