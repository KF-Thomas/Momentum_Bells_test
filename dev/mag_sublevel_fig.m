%energy level diagram

levels = {'2^3S_1','m_J=-1','m_J=0','m_J=+1','2^3P_0','2^1P_1','3^3S_1','2^3P'};
    
energy = [20.61577,20.61577-0.5,20.61577,20.61577+0.5];
energy_adj = [20.61577,20.61577-0.5,20.61577,20.61577+0.5];
state_pos = [1,1.9,1.9,1.9];
lin_width = 0.25;
y_shift = 0.1;
x_shift = 0.15;
font_size=18.9;
stfig('energy diagram')
clf
for ii = 1:length(state_pos)
    level_str = levels{ii};
    pos = state_pos(ii);
    energy_current = energy_adj(ii);
    plot(pos+[-lin_width,lin_width],energy_current.*[1,1],'k-','LineWidth',2.0)
    if ii<7 && ii>3
        text(pos-x_shift,energy_current-y_shift,['\(',level_str,'\)'],'interpreter','latex','FontSize',font_size)
    else
        text(pos-x_shift,energy_current-y_shift,['\(',level_str,'\)'],'interpreter','latex','FontSize',font_size)
    end
    hold on
end
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
%add the arrows
% quiver(state_pos(2),energy_adj(2),state_pos(8)-state_pos(2),energy_adj(8)-energy_adj(2),1,'filled','color',[0.4940 0.1840 0.5560])
% quiver(state_pos(2),energy_adj(2),state_pos(9)-state_pos(2),energy_adj(9)-energy_adj(2),0,'filled')
% quiver(state_pos(9),energy_adj(9),state_pos(8)-state_pos(9),energy_adj(8)-energy_adj(9),0,'filled')
% quiver(state_pos(1),energy_adj(1),state_pos(2)-state_pos(1),energy_adj(2)-energy_adj(1),0,'filled')

plot([state_pos(1)+lin_width,state_pos(2)-lin_width],[energy_adj(1),energy_adj(2)],'k-.')
plot([state_pos(1)+lin_width,state_pos(3)-lin_width],[energy_adj(1),energy_adj(3)],'k-.')
plot([state_pos(1)+lin_width,state_pos(4)-lin_width],[energy_adj(1),energy_adj(4)],'k-.')
box off

set(gca,'Visible','off')
%scatter(state_pos,energy_adj,'kx')