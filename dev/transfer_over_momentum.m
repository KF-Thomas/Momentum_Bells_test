num_pts=100;
opts.control=0;
opts.timer=0;
%gaussian pulse
% [17 0.0235 0.27 -0.085*2*2*pi, 0.00]
% [17 0.03 0.31 -0.0849*2*2*pi, 0.00] %cost = 4.299
% [17 0.0245 0.29 -0.0849*2*2*pi, 0.00]
% [17 0.02 0.25 -0.0849*2*2*pi, 0.00] %cost = 4.972
% [17 0.05 0.405 -0.0849*2*2*pi, 0.00] %cost = 3.452

%sinc squared pulse
% [17 45 0.325 -0.0849*2*2*pi, 0.00] %cost = 5.226
% [17 10 0.325 -0.0849*2*2*pi, 0.00] %cost = 3.821

%sinc pulse
% [17 6 0.58 -0.0849*2*2*pi, 0.00] %cost = 1.531
% [17 6 0.508 -0.0849*2*2*pi, 0.00] %cost = 1.297
% [17.5 5.8 0.508 -0.0849*2*2*pi, 0.00] %cost = 1.257
% [17.5 5.3 0.59 -0.0849*2*2*pi, 0.00] %cost = 1.245

%rectangle pulse
% [10 6 0.54 -0.0849*2*2*pi, 0.00] %cost = 2.945
% [17 6 0.58 -0.0849*2*2*pi, 0.00] %cost = 3.143

%sinc times gaussian pulse

opts.time_span = 35;
ys=transfer_percentage([17.5 4 0.6 -0.0849*2*2*pi, 0.00, 0.025],[9,10,11],opts,linspace(0,4,num_pts).*-4.101078618245948e+06);
%% cost
k = linspace(0,4,num_pts).*-4.101078618245948e+06;
mask_t =  logical(k<0.25.*-4.101078618245948e+06) & ...
        logical(k>0.75.*-4.101078618245948e+06);
    mask_b =  (k<2.25.*-4.101078618245948e+06) & ...
        (k>3.75.*-4.101078618245948e+06);
    cost = trapz(abs(k(mask_t))./1e6,1-abs(ys(mask_t,1)).^2)+...
        trapz(abs(k(mask_b))./1e6,1-abs(ys(mask_b,3)).^2)


%%
stfig('transfer over momentum space')
plot(linspace(0,4,num_pts),abs(ys).^2,'linewidth',2)
xlabel('$-k/k_0$')
ylabel('transfer percentage')
legend('n=-1','n=0','n=+1')
ylim([0 1])