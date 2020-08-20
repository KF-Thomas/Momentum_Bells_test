this_folder = fileparts(which(mfilename));
idcs   = strfind(this_folder,'\');
this_folder = this_folder(1:idcs(end)-1);
addpath(genpath(this_folder));
core_folder = fullfile(fileparts(this_folder), 'Core_BEC_Analysis\');
addpath(genpath(core_folder));
set(groot, 'DefaultTextInterpreter', 'latex')

cost_opts.diffraction_order = 9;
cost_opts.orders = [-1,0]; %the difraction orders we want top optimize transfer to

cost_fun = @(b) Raman_Nath_Cost(b,cost_opts);
% b0 = [2.0546 0.0595 1.1221 0.1682];%[1,2,1.8,-0.1]; %start of gaussian pulse, envoloping, amplitude, detuning
% b0 = [0.4572 0.0281 1.0723 0.1674];
% b0 = [4.4794 0.0122 1.2537 0.2300];
% b0 = [5 0.0138 0.7983 -0.8063];
% b0 = [5.0000    0.1000    0.2918   -0.4707];
b0 = [4.9999    0.1000    0.2887   -0.5060];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1, 0.1, 0.2, -0.9];
ub = [5, 6, 5, 0.9];
nonlcon = [];
% options_fminsearch = optimset('PlotFcns',@optimplotfval);
% options = optimoptions('PlotFcns',@optimplotfval,'patternsearch','Cache','on')
options_fmincon = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',...
    1500,'FiniteDifferenceType','forward');


% b = fminsearch(cost_fun,b0,options_fminsearch);
b_opt = fmincon(cost_fun,b0,A,b,Aeq,beq,lb,ub,nonlcon,options_fmincon);%,options);
% b_opt = patternsearch(cost_fun,b0,A,b,Aeq,beq,lb,ub,nonlcon);%,options);

[costmin cost] = Raman_Nath_Cost(b_opt,cost_opts);
costmin