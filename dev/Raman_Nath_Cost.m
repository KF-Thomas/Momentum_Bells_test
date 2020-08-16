function [final_cost, cost] = Raman_Nath_Cost(b,cost_opts)
solver_opts.plots = false;
solver_opts.timer = false;
solver_opts.diffraction_order = cost_opts.diffraction_order;
[t,y] = Raman_Nath_Solver(b,solver_opts);
indxs = cost_opts.orders+cost_opts.diffraction_order+1;
nmodes = length(indxs);
Norm = sum(abs(y).^2,2);
cost = 0;
for ii = indxs
    cost = cost + abs(abs(y(:,ii)).^2-Norm/nmodes);
end
final_cost = cost(end);
end