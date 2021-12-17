%%
% load('C:\Users\BEC Machine\cloudstor\PROJECTS\Momentum_Bells_test\sims\20210803_k=0,-1,-2_rt_scan_mid_trap_equal_delay\g2amp_vs_dk.mat')
load('C:\Users\BEC Machine\cloudstor\PROJECTS\Momentum_Bells_test\sims\g2amp_vs_dk_20210824')
hebec_constants
indx_z = length(k_z_12{1});
indx_vec = (1:indx_z)';
dk_vec_z = [1e-3.*indx_vec,6e-3.*ones(indx_z,1),2.4e-2.*ones(indx_z,1)];

indx_x = length(k_x_12{1});
indx_vec = (1:indx_x)';
dk_vec_x = [1.68e-2.*ones(indx_x,1),1e-3.*indx_vec,2.4e-2.*ones(indx_x,1)];

indx_y = length(k_y_12{1});
indx_vec = (1:indx_y)';
dk_vec_y = [1.68e-2.*ones(indx_y,1),6e-3.*ones(indx_y,1),1e-3.*indx_vec];

k_min = 5*[1,1,1].*1e-3;%[4.2e-3,1.5e-3,6e-3].*4.*0.05;
k_max = [4.2e-3,1.5e-3,6e-3].*4.*5;
z_mask = dk_vec_z(:,1)>k_min(1);
x_mask = dk_vec_x(:,2)>k_min(2);
y_mask = dk_vec_y(:,3)>k_min(3);

dk_vec = [dk_vec_z(z_mask,:);dk_vec_x(x_mask,:);dk_vec_y(y_mask,:)];

colors_main=[[220,80,80];[60,60,219]./1.3;[219,60,219]./1.5;[60,220,135]./1.5];
colors_main=colors_main./255;

colors_main_2=[[210,113,80];[60,220,135]./1.5;[88,113,219]./1.5;[60,60,219]./1.5]./1.5;
colors_main_2=colors_main_2./255;
%%
jj = 3;
E_z = (k_z_12{jj}(z_mask)+k_z_34{jj}(z_mask)-k_z_14{jj}(z_mask)-k_z_23{jj}(z_mask))./...
    (k_z_12{jj}(z_mask)+k_z_34{jj}(z_mask)+k_z_14{jj}(z_mask)+k_z_23{jj}(z_mask));
E_x = (k_x_12{jj}(x_mask)+k_x_34{jj}(x_mask)-k_x_14{jj}(x_mask)-k_x_23{jj}(x_mask))./...
    (k_x_12{jj}(x_mask)+k_x_34{jj}(x_mask)+k_x_14{jj}(x_mask)+k_x_23{jj}(x_mask));
E_y = (k_y_12{jj}(y_mask)+k_y_34{jj}(y_mask)-k_y_14{jj}(y_mask)-k_y_23{jj}(y_mask))./...
    (k_y_12{jj}(y_mask)+k_y_34{jj}(y_mask)+k_y_14{jj}(y_mask)+k_y_23{jj}(y_mask));

mdl_full_z = @(b,x) rob_E([x(:,1)./(b(2)),x(:,2)./(4.63796173807889e-03),x(:,3)./(1.45518272470288e-02)]',[b(1)],0)';
mdl_fit_g2_z=fitnlm(dk_vec_z(z_mask,:),E_z',mdl_full,[10,[4.2e-3.*4]])

mdl_full_x = @(b,x) rob_E([x(:,1)./(7.30336375932467e-03),x(:,2)./(b(2)),x(:,3)./(1.45518272470288e-02)]',[b(1)],0)';
mdl_fit_g2_x=fitnlm(dk_vec_x(x_mask,:),E_x',mdl_full,[10,[4.2e-3.*4]])

mdl_full_y = @(b,x) rob_E([x(:,1)./(7.30336375932467e-03),x(:,2)./(4.63796173807889e-03),x(:,3)./(b(2))]',[b(1)],0)';
mdl_fit_g2_y=fitnlm(dk_vec_y(y_mask,:),E_y',mdl_full,[10,[4.2e-3.*4]])

%%
jj=3; %select which phase we want to look at
port_list = {'14','23','12','34'};
styl = {'r','b','m','g'};
do_ports = 1:4;%[1,2,4];%
for mm = 1:length(do_ports)
    kk = do_ports(mm);
current_data_z = eval(['k_z_',port_list{kk}]);
current_data_x = eval(['k_x_',port_list{kk}]);
current_data_y = eval(['k_y_',port_list{kk}]);
wz = eval(['w_z_',port_list{kk}]);
wx = eval(['w_x_',port_list{kk}]);
wy = eval(['w_y_',port_list{kk}]);

w_vec = [wz{jj}(z_mask);wx{jj}(x_mask);wy{jj}(y_mask)];
y = [current_data_z{jj}(z_mask);current_data_x{jj}(x_mask);current_data_y{jj}(y_mask)];

mdl_full = @(b,x) rob_g([x(:,1)./(b(2)),x(:,2)./(b(3)),x(:,3)./(b(4))]',[b(1);1;b(5)],b(6))';
mdl_fit_g2_full=fitnlm(dk_vec,y,mdl_full,[0,[4.2e-3,1.5e-3,6e-3].*4,0,0],'Weight',1./w_vec.^0.5)

stfig('g2 vs dk');
if mm == 1
clf
end
xp=linspace(0,5e-2,1000)';
subplot(3,1,1)
kp=[xp,6e-3.*ones(length(xp),1),2.4e-2.*ones(length(xp),1)];

[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_full,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
% [ysamp_val_2,ysamp_ci_2]=predict(mdl_fit_g2_fixed,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));
%errorbar(lambda,top_corr_bb_vec(:,1),wt(:,1),'kx')

hg(kk)=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

% hfixed=plot(xp,ysamp_val_2,'r','LineWidth',1.5);
% plot(xp,ysamp_ci_2,'color',[1,1,1].*0.5)

errorbar(dk_vec_z(z_mask,1).*const.mhe./const.hb.*1e-6,current_data_z{jj}(z_mask),wz{jj}(z_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_z(z_mask,1).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_z$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])
subplot(3,1,2)
kp=[1.68e-2.*ones(length(xp),1),xp,2.4e-2.*ones(length(xp),1)];
[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_full,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
hfree=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

errorbar(dk_vec_x(x_mask,2).*const.mhe./const.hb.*1e-6,current_data_x{jj}(x_mask),wx{jj}(x_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_x(x_mask,2).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_x$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])

subplot(3,1,3)
kp=[1.68e-2.*ones(length(xp),1),6e-3.*ones(length(xp),1),xp];
[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_full,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
hfree=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

errorbar(dk_vec_y(y_mask,3).*const.mhe./const.hb.*1e-6,current_data_y{jj}(y_mask),wy{jj}(y_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_y(y_mask,3).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_y$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])
end
legend([hg(1) hg(2) hg(3) hg(4)],{'$g^{(2)}_{14}$','$g^{(2)}_{23}$','$g^{(2)}_{12}$', '$g^{(2)}_{34}$'})

%%
jj=1;
w_all = [];
y_all = [];
phi_offset = [];
dk_all = [];
for kk = 1:4
current_data_z = eval(['k_z_',port_list{kk}]);
current_data_x = eval(['k_x_',port_list{kk}]);
current_data_y = eval(['k_y_',port_list{kk}]);
wz = eval(['w_z_',port_list{kk}]);
wx = eval(['w_x_',port_list{kk}]);
wy = eval(['w_y_',port_list{kk}]);

w_vec = [wz{jj}(z_mask);wx{jj}(x_mask);wy{jj}(y_mask)];
y = [current_data_z{jj}(z_mask);current_data_x{jj}(x_mask);current_data_y{jj}(y_mask)];

if kk>2
    flip = 0;
else
    flip = 1;
end
phi_currrent = flip.*pi.*ones(size(w_vec));

w_all = [w_all;w_vec];
y_all = [y_all;y];
phi_offset = [phi_offset;phi_currrent];
dk_all = [dk_all;dk_vec];
end
mdl_full = @(b,x) rob_g([x(:,1)./(b(2)),x(:,2)./(b(3)),x(:,3)./(b(4))]',[b(1);1;b(5)],b(6)+x(:,4)')';
mdl_fit_g2_all=fitnlm([dk_all,phi_offset],y_all,mdl_full,[20,[4.2e-3,1.5e-3,6e-3].*4,0,0],'Weight',1./w_all.^0.5)
% 
% mdl_full = @(b,x) rob_g([x(:,1)./(b(2)),x(:,2)./(b(3)),x(:,3)./(b(4))]',[b(1);1;0],b(5)+x(:,4)')';
% mdl_fit_g2_all=fitnlm([dk_all,phi_offset],y_all,mdl_full,[20,[4.2e-3,1.5e-3,6e-3].*4,0],'Weight',1./w_all.^0.5)
% % 
%%
stfig('g2 vs dk full fit');
for kk = 1:4
    current_data_z = eval(['k_z_',port_list{kk}]);
current_data_x = eval(['k_x_',port_list{kk}]);
current_data_y = eval(['k_y_',port_list{kk}]);
wz = eval(['w_z_',port_list{kk}]);
wx = eval(['w_x_',port_list{kk}]);
wy = eval(['w_y_',port_list{kk}]);

w_vec = [wz{jj}(z_mask);wx{jj}(x_mask);wy{jj}(y_mask)];
y = [current_data_z{jj}(z_mask);current_data_x{jj}(x_mask);current_data_y{jj}(y_mask)];

if kk == 1
clf
end
if kk>2
    flip = 0;
else
    flip = 1;
end
xp=linspace(0,5e-2,1000)';
subplot(3,1,1)
kp=[xp,6e-3.*ones(length(xp),1),2.4e-2.*ones(length(xp),1),flip.*ones(length(xp),1)];

[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_all,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
% [ysamp_val_2,ysamp_ci_2]=predict(mdl_fit_g2_fixed,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2)));
%errorbar(lambda,top_corr_bb_vec(:,1),wt(:,1),'kx')

hg2(kk)=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

% hfixed=plot(xp,ysamp_val_2,'r','LineWidth',1.5);
% plot(xp,ysamp_ci_2,'color',[1,1,1].*0.5)

errorbar(dk_vec_z(z_mask,1).*const.mhe./const.hb.*1e-6,current_data_z{jj}(z_mask),wz{jj}(z_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_z(z_mask,1).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_z$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])
subplot(3,1,2)
kp=[1.68e-2.*ones(length(xp),1),xp,2.4e-2.*ones(length(xp),1),flip.*ones(length(xp),1)];
[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_all,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
hfree=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

errorbar(dk_vec_x(x_mask,2).*const.mhe./const.hb.*1e-6,current_data_x{jj}(x_mask),wx{jj}(x_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_x(x_mask,2).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_x$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])

subplot(3,1,3)
kp=[1.68e-2.*ones(length(xp),1),6e-3.*ones(length(xp),1),xp,flip.*ones(length(xp),1)];
[ysamp_val,ysamp_ci]=predict(mdl_fit_g2_all,kp,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
hfree=plot(xp.*const.mhe./const.hb.*1e-6,ysamp_val,styl{kk},'LineWidth',1.5);
hold on
% plot(xp.*const.mhe./const.hb.*1e-6,ysamp_ci,'color',[1,1,1].*0.5)

errorbar(dk_vec_y(y_mask,3).*const.mhe./const.hb.*1e-6,current_data_y{jj}(y_mask),wy{jj}(y_mask).^0.5,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(kk,:),...
    'MarkerFaceColor',colors_main_2(kk,:),'LineWidth',2.5)
xlim([0 1.05.*max(dk_vec_y(y_mask,3).*const.mhe./const.hb.*1e-6)])
ylabel('$g^{(2)}_{bb}(0)$','interpreter','latex')
xlabel('$\Delta k_y$ ($\mu$m$^{-1}$)')
% legend([hfree hfixed],'Free','Fixed')
set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
font_size_global = 15;
set(gca,'fontsize',font_size_global)
ylim([0 10])
end