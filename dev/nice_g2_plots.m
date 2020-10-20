function nice_g2_plots(corr_opts,ports)
port_labels = {'g14','g12','g23','g34'};
centers = 'rad_centers';

colors_main=[[88,113,219];[60,220,180]./1.75;[88,113,219]./1.7]; %[88,113,219]%[96,144,201]
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=13;

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);
%which model to use

for ii = 1:4
port = port_labels{ii};
corrs = ports.(port).norm_g2;
stfig([port,' distribution']);
clf
    if corr_opts.calc_err
        errorbar(corrs.(centers).*1e3,corrs.g2_amp,xg2_err,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_main(2,:),'LineWidth',2.5)
    else
        plot(corrs.(centers).*1e3,corrs.g2_amp,'o','MarkerSize',5,'MarkerFaceColor',colors_detail(1,:),'MarkerEdgeColor',colors_main(2,:))
    end
    ylabel('$g^{(2)}(\Delta v)$','interpreter','latex')
    xlabel('$\Delta v$ (mm/s)','interpreter','latex')
    if corr_opts.fit
        hold on
        xx = linspace(min([max(corrs.(centers)),0]),max(corrs.(centers)),3e3)';
        if ~corr_opts.calc_err
            fit =ports.(port).fit;
            [ypred,ypredci] = predict(fit,xx,'Simultaneous',true);
            
            curve1 = ypredci(:,1)';
            curve2 = ypredci(:,2)';
            x1 = (xx.*1e3)';
            x2 = [x1, fliplr(x1)];
            inBetween = [curve1, fliplr(curve2)];
            h = fill(x2, inBetween, 'g');
            h.FaceColor = [0.31 0.31 0.32].*2;
            h.FaceAlpha = 0.5;
            
            plot(xx.*1e3,ypred,'b-', xx.*1e3,ypredci,'r-');
        else
            plot(xx.*1e3,fun1d(b,xx),'b-')
            plot(xx.*1e3,fun1d(b-b_unc,xx),'r-')
            plot(xx.*1e3,fun1d(b+b_unc,xx),'r-')
        end
    end
    set(gca,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.1)
ax = gca;
ax.XAxis.TickLabelFormat= '\\textbf{%g}';
ax.YAxis.TickLabelFormat= '\\textbf{%g}';
set(gca,'fontsize',font_size_global)
xlim([0,max(corrs.(centers).*1e3)])
end