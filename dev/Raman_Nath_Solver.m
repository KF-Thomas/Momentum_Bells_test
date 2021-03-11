function [t,y] = Raman_Nath_Solver(b,opts,k)
if nargin<3
    k=0;
end
if opts.timer
    tic
end
%% USER INPUTS
if isfield(opts,'diffraction_order')
    diffraction_order = opts.diffraction_order;
else
    diffraction_order = 9;% number of orders diffraction (either side of zero) we wish to consider nth ifrraction order + 1
end
if isfield(opts,'time_span')
    tspan = [0 opts.time_span];
else
    tspan = [0 35];% time span in mus
end

%% SYSTEM CONSTANTS
hebec_constants
laser_wavelength = 1083.33e-9; %wavelength of the laser in m

%% GENERATED PARAMETERS
N = diffraction_order+1;
recoil_freqency = const.hb*(sqrt(2)*2*pi/laser_wavelength)^2/(2*pi*2*const.mhe);% the frequency of the photon momentum
wr = recoil_freqency*2*pi/4/1e6;% recoil frequency diveded by 4 in rad*MHz
% c0 = [zeros(N-1,1);-1i;zeros(N-1,1)];
c0 = [zeros(N-2,1);0.5;0.5;zeros(N-1,1)];
%% SOLVE THE RAMAN-NATH EQUATIONS
func = @(t,c) Raman_Nath(t,c,N,wr,b,k);
[t,y] = ode45(func,tspan,c0);
if opts.timer
    toc
end

%% PLOT OUT RESULTS
if isfield(opts,'plot') && opts.plot
    stfig('Raman-Nath simulation');
    clf
    plot(t,abs(y(:,(-2:2)+N)).^2,'Linewidth',2.2)
    xlabel('time ($\mu$s)','interpreter','latex')
    ylabel('$|C_n|^2$','interpreter','latex')
    test=(-2:2);
    legend(string(test))
end

end



