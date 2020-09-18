num_pts=200;
opts.control=0;
opts.timer=0;
%gaussian pulse
% [17 0.0235 0.27 -0.085*2*2*pi, 0.00]
% [17 0.03 0.31 -0.0849*2*2*pi, 0.00] %cost = 4.299
% [17 0.0245 0.29 -0.0849*2*2*pi, 0.00]
% [17 0.02 0.25 -0.0849*2*2*pi, 0.00] %cost = 4.972
% [17 0.05 0.405 -0.0849*2*2*pi, 0.00] %cost = 3.452

%super gaussian pulse
% [17 0.05^2 0.39 4 -0.0849*2*2*pi 0.000] %cost = 3.0254
% [17 0.05^2 0.375 4 -0.0849*2*2*pi 0.000] %cost = 3.0107
% [17 0.06^2 0.375 4 -0.0849*2*2*pi 0.000] %cost = 2.914
% [17 0.06^2 0.416 4 -0.0849*2*2*pi 0.000] %cost = 2.838
% [17 0.082^3 0.455 6 -0.0849*2*2*pi 0.000] %cost = 2.738
% [17 0.1 0.5 2 -0.0849*2*2*pi 0.000]

%sinc squared pulse
% [17 45 0.325 -0.0849*2*2*pi, 0.00] %cost = 5.226
% [17 10 0.325 -0.0849*2*2*pi, 0.00] %cost = 3.821

%sinc pulse
% [17 6 0.58 -0.0849*2*2*pi, 0.00] %cost = 1.531
% [17 6 0.508 -0.0849*2*2*pi, 0.00] %cost = 1.297
% [17.5 5.8 0.508 -0.0849*2*2*pi, 0.00] %cost = 1.257
% [17.5 5.3 0.59 -0.0849*2*2*pi, 0.00] %cost = 1.245

% [17.5 6.0 0.5151 -0.0849*2*2*pi 0.00]

%rectangle pulse
% [10 6 0.54 -0.0849*2*2*pi, 0.00] %cost = 2.945
% [17 6 0.58 -0.0849*2*2*pi, 0.00] %cost = 3.143

%hamming sinc pulse
% [17.5 6.0 0.5151 -0.0849*2*2*pi 0.00]

opts.time_span = 35;
ks = [linspace(0,2,num_pts).*-4.101078618245948e+06 linspace(2,4,num_pts).*-4.101078618245948e+06];
ys=transfer_percentage([17 0.05 0.405 -0.0849*2*2*pi, 0.00],[9,10,11],opts,ks);
%% cost
k = linspace(0,4,num_pts).*-4.101078618245948e+06;
mask_t =  logical(k<0.25.*-4.101078618245948e+06) & ...
        logical(k>0.75.*-4.101078618245948e+06);
    mask_b =  (k<2.25.*-4.101078618245948e+06) & ...
        (k>3.75.*-4.101078618245948e+06);
    cost = trapz(abs(k(mask_t))./1e6,1-abs(ys(mask_t,1)).^2)+...
        trapz(abs(k(mask_b))./1e6,1-abs(ys(mask_b,3)).^2)


%%
num_pts=400;
stfig('transfer over momentum space')
clf
% hold on
plot(ks./-4.101078618245948e+06,abs(ys(:,:)).^2,'linewidth',2)
% hold on

xlabel('$-k/k_0$')
ylabel('transfer percentage')
legend('n=-1','n=0','n=+1')
ylim([0 1])
grid on
stfig('phase over momentum space')
clf
plot(linspace(0,4,num_pts),unwrap((angle(ys(:,[1,3])))+(2*pi)/2.*(1-sign(angle(ys(:,[1,3])))))./(2*pi),'linewidth',2)
% hold off

xlabel('$-k/k_0$')
ylabel('phase/2$\pi$')
legend('n=-1','n=+1')
% ylim([0 1])
grid on

stfig('phase differnce over momentum')
clf
kt=linspace(0,2,num_pts/2);
% dif_1 = angle(ys(1:(num_pts/2),1))-angle(ys((num_pts/2):end-1,3));
dif_2 = angle(ys(1:(num_pts/2),1))-angle(ys((num_pts/2+1):end,3));
% plot(kt,mod(dif_1,2*pi)./(2*pi),'linewidth',2)
hold on
plot(kt,unwrap(mod(dif_2,2*pi))./(2*pi),'linewidth',2)
xlabel('$-k/k_0$')
ylabel('$\delta$phase/2$\pi$')
%%
opts.control=1;
num_pts=220;
amp_vec=linspace(0.0,4.5,num_pts);
% linspace(0,4,num_pts).*
% [17 0.1 0.5 2 -0.0849*2*2*pi 0.000]
% b_matrix = [17.*ones(num_pts,1) 0.07.*ones(num_pts,1) linspace(0.0,1.4,num_pts)' 2.*ones(num_pts,1) -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)];
b_matrix_gauss = [17.*ones(num_pts,1) 0.05.*ones(num_pts,1) amp_vec' -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)];
% b_matrix_sinc = [17.5.*ones(num_pts,1) 6.*ones(num_pts,1) amp_vec' -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)]; %cost = 1.245

ybs=transfer_percentage(b_matrix_gauss,[9,10,11],opts,-4.101078618245948e+06);
%
stfig('transfer over parameter space')
clf
% hold on
plot(amp_vec,abs(ybs).^2,'linewidth',2)
grid on
xlabel('amp')
ylabel('transfer percentage')
legend('n=-1','n=0','n=+1')
ylim([0 1])
%%
opts.control=1;
alpha_vec = linspace(0.01,1,num_pts)';
% linspace(0,4,num_pts).*
% [17 0.1 0.5 2 -0.0849*2*2*pi 0.000]
% b_matrix = [17.*ones(num_pts,1) 0.07.*ones(num_pts,1) linspace(0.0,1.4,num_pts)' 2.*ones(num_pts,1) -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)];
% b_matrix = [17.5.*ones(num_pts,1) alpha_vec 0.5151.*ones(num_pts,1) -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)];
b_matrix_gauss = [17.0.*ones(num_pts,1) alpha_vec 0.411.*ones(num_pts,1) -0.0849*2*2*pi.*ones(num_pts,1) 0.00.*ones(num_pts,1)];

yAs=transfer_percentage(b_matrix_gauss,[9,10,11],opts,-4.101078618245948e+06);
%
stfig('transfer over parameter space (alpha)')
clf
% hold on
plot(alpha_vec,abs(yAs).^2,'linewidth',2)

xlabel('alpha')
ylabel('transfer percentage')
legend('n=-1','n=0','n=+1')
ylim([0 1])