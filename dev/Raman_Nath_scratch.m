opts.plot=1;
opts.control=0;
opts.timer=0;
opts.time_span = 50;
params = [25 sqrt(1/(2*0.05)) 0.405 0, 0.00];
%[17 0.0235 0.27 -0.085*2*2*pi, 0.0]
[t,y] = Raman_Nath_Solver(params,opts,-4.101078618245948e+06);
%%
stfig('Raman-Nath phase')
    clf
% hold on
    plot(t,angle(y(:,[11,9])),'Linewidth',2.2)
    xlabel('time ($\mu$s)','interpreter','latex')
    ylabel('$ang(C_n)$','interpreter','latex')
    
        stfig('Raman-Nath mag');
    clf
    plot(t,abs(y(:,[11,9])).^2,'Linewidth',2.2)
    xlabel('time ($\mu$s)','interpreter','latex')
    ylabel('$|C_n|^2$','interpreter','latex')
    %%
    y_m1 = y(:,9);
    y_1 = y(:,11);
t_m1=t;
figure(111)
clf
hold on
plot(t,mod(angle(y_1),2*pi))
plot(t_m1,mod(angle(y_m1),2*pi))