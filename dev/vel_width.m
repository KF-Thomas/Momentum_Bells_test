% momentum width scratch

omega=10; %in hz*(2*pi)
tfall=0.416;
widthff=@(omega) sqrt( ...
                        (sqrt(const.hb./(2*omega*(2*pi).*const.mhe))).^2  ...
                        +(tfall*sqrt(omega*(2*pi)*const.hb./(2*const.mhe)  )   ).^2 ...
                        );
% width_x=@(omega) sqrt(const.hb./(2*omega*(2*pi).*const.mhe))     ;               
%                     
% widthff=@(omega) sqrt( ...
%     (width_x(omega)).^2  ...
%     +(tfall*(const.hb/2)*(1/const.mhe)*(1./width_x(omega))   ).^2 ...
%     );


[omega_min,width_min]=fminsearch(widthff,0.6)

xsamp=linspace(0.1,2,1e3);
stfig('detector width')
plot(xsamp,widthff(xsamp)*1e6)
xlabel('trap freq (Hz)')
ylabel('detector width ($\mu$m)')



%%
widthff(5.8)*1e6

fminsearch(@(x) abs(widthff(x)-120e-6),3)


%%
combined_res=2*120e-6
%combined_res=sqrt(2)*120e-6
[omega,width]=fminsearch(@(x) abs(sqrt(widthff(x)^2+120e-6^2)-combined_res),1)

