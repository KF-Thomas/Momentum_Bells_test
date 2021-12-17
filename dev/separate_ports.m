function [left_port, right_port] = separate_ports(data,opts_ports)
%takes in the halo in velocity space and returns the two separate ports
%indexes [z , x , y] = [3 , 1 , 2]
% split along the x-axis
left_port.unnorm = cell(size(data.counts_vel));
right_port.unnorm = cell(size(data.counts_vel));

if ~isfield(opts_ports,'norm') || opts_ports.norm
    left_port.norm = cell(size(data.counts_vel));
    right_port.norm = cell(size(data.counts_vel));
end
for ii = 1:size(data.counts_vel,1)
    current_shot = data.counts_vel{ii};
    mask = current_shot(:,2)<0; % port mask
    if isfield('opts_ports','pol_lims')
    mask_pol = (atan2(current_shot(:,2),current_shot(:,3))<opts_ports.pol_lims(2) &...
        atan2(current_shot(:,2),current_shot(:,3))>opts_ports.pol_lims(1)) |...
        (atan2(current_shot(:,2),current_shot(:,3))<(opts_ports.pol_lims(2)-pi) &...
        atan2(current_shot(:,2),current_shot(:,3))>(opts_ports.pol_lims(1))-pi);% polar mask
    mask = mask & mask_pol;
    end
    mask_norm = current_shot(:,2)<0;
    left_port.unnorm{ii} = current_shot(mask,:);
    right_port.unnorm{ii} = current_shot(~mask,:);
    if ~isfield(opts_ports,'norm') || opts_ports.norm
        current_shot_norm = data.counts_vel_norm{ii};
        left_port.norm{ii} = current_shot_norm(mask_norm,:);
        right_port.norm{ii} = current_shot_norm(~mask_norm,:);
    end
end

end