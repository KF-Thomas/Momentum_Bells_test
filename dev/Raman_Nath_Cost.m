function [final_cost, cost] = Raman_Nath_Cost(b,cost_opts)
solver_opts.plots = false;
solver_opts.timer = false;
solver_opts.diffraction_order = cost_opts.diffraction_order;

indxs = cost_opts.orders+cost_opts.diffraction_order+1;
nmodes = length(indxs);

if strcmp(cost_opts.goal,'transfer')
    [t,y] = Raman_Nath_Solver(b,solver_opts);
    Norm = sum(abs(y).^2,2);
    cost = 0;
    for ii = indxs
        cost = cost + abs(abs(y(:,ii)).^2-Norm/nmodes);
    end
elseif strcmp(cost_opts.goal,'mirror')
    solver_opts.control = 0;
    num_pts = 100;
    ords = cost_opts.diffraction_order:(cost_opts.diffraction_order+2);
    k = linspace(0,4,num_pts).*-4.101078618245948e+06;
    ys=transfer_percentage(b,ords,solver_opts,k);
    mask_t =  k<0.25.*-4.101078618245948e+06 && ...
        k>0.75.*-4.101078618245948e+06;
    mask_b =  k<2.25.*-4.101078618245948e+06 && ...
        k>3.75.*-4.101078618245948e+06;
    cost = sum(1-abs(ys(mask_t,cost_opts.diffraction_order)).^2)+...
        sum(1-abs(ys(mask_b,cost_opts.diffraction_order+2)).^2);
end
final_cost = cost(end);
end