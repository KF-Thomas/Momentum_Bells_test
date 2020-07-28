function [left_port, right_port] = separate_ports(data,opts_ports)
%takes in the halo in velocity space and returns the two separate ports
%indexes [z , x , y] = [3 , 1 , 2]
% split along the x-axis
left_port = cell(size(data.counts_vel));
right_port = cell(size(data.counts_vel));
for ii = 1:size(data.counts_vel,1)
    current_shot = data.counts_vel{ii};
    mask = current_shot(:,2)<0;
    left_port{ii} = current_shot(mask,:);
    right_port{ii} = current_shot(~mask,:);
end

end