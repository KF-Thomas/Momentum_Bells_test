[X,Y,Z] = sphere(100);
colormap(jet)
R = 1;
hold on
sn = surf(R*X,R*Y, R*(Z+1),0.*X+1);
ss = surf(R*X,R*Y,R*(Z-1),0.*X-1);
sn.LineStyle = 'none';
ss.LineStyle = 'none';
sn.FaceAlpha = 0.3;
ss.FaceAlpha = 0.3;
view([30 40])
axis equal
%% 
h= figure;
h.Color = [.2 .2 .2];
ax = axes;
ax.Color = [.2 .2 .2];
ax.XColor = [1 1 1];
ax.YColor = [1 1 1];
ax.ZColor = [1 1 1];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
ax.GridColor = [1 1 1];

set(gca,'nextplot','replacechildren');
loops = 60;
R = @(j) j/(3*loops);
M(loops) = struct('cdata',[],'colormap',[]);
%h.Visible = 'off'
v = VideoWriter('GrowthSphere.avi');
v.FrameRate = 12;
open(v);
for j = 1:loops
    points = [R(j)/2 R(j)/2 R(j)*(1+1/sqrt(2)); -R(j)/2 -R(j)/2 R(j)*(1-1/sqrt(2)); R(j)/2 R(j)/2 R(j)*(-1+1/sqrt(2)); -R(j)/2 -R(j)/2 -R(j)*(1+1/sqrt(2))]; 
    hold off
    sn = surf(ax,X*R(j),Y*R(j), R(j)*(Z+1),0.*X+1);
    hold on
    ss = surf(ax,R(j)*X,R(j)*Y,R(j)*(Z-1),0.*X-1);
    scatter3(ax,points(:,1),points(:,2),points(:,3),10*R(j),'x')
    plot3(ax,points([1 2 4 3 1],1),points([1 2 4 3 1],2),points([1 2 4 3 1],3),'-.k')
    ylim([-1 1])
    xlim([-1 1])
    zlim([-2 2])
    sn.LineStyle = 'none';
    ss.LineStyle = 'none';
    sn.FaceAlpha = 0.3;
    ss.FaceAlpha = 0.3;
    
    ax.Color = [.2 .2 .2];
    ax.XColor = [1 1 1];
    ax.YColor = [1 1 1];
    ax.ZColor = [1 1 1];
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    ax.GridColor = [1 1 1];
    view([30 40])
    daspect([1 1 1]);
    drawnow
    M(j) = getframe(gcf);
    writeVideo(v,M(j));
end
for j=loops+1:2*loops
    points = [R(j)/2 R(j)/2 R(j)*(1+1/sqrt(2))-(2/3)*(j-loops)/loops; -R(j)/2 -R(j)/2 R(j)*(1-1/sqrt(2))-(2/3)*(j-loops)/loops; R(j)/2 R(j)/2 R(j)*(-1+1/sqrt(2))+(2/3)*(j-loops)/loops; -R(j)/2 -R(j)/2 -R(j)*(1+1/sqrt(2))+(2/3)*(j-loops)/loops]; 
    hold off
    axes(ax)
    sn = surf(X*R(j),Y*R(j), R(j)*(Z+1)-(2/3)*(j-loops)/loops,0.*X+1);
    hold on
    ss = surf(R(j)*X,R(j)*Y,R(j)*(Z-1)+(2/3)*(j-loops)/loops,0.*X-1);
    scatter3(points(:,1),points(:,2),points(:,3),10*R(j),'x')
    plot3(points([1 2 4 3 1],1),points([1 2 4 3 1],2),points([1 2 4 3 1],3),'-.k')
    ylim([-1 1])
    xlim([-1 1])
    zlim([-2 2])
    sn.LineStyle = 'none';
    ss.LineStyle = 'none';
    sn.FaceAlpha = 0.3;
    ss.FaceAlpha = 0.3;
    view([30 40])
    daspect([1 1 1]);
    drawnow
    M(j) = getframe(gcf);
    writeVideo(v,M(j));
end
for j=2*loops+1:3*loops
    points = [R(j)/2 R(j)/2 R(j)*(1+1/sqrt(2))-(2/3)*(2*j-3*loops)/loops; -R(j)/2 -R(j)/2 R(j)*(1-1/sqrt(2))-(2/3)*(2*j-3*loops)/loops; R(j)/2 R(j)/2 R(j)*(-1+1/sqrt(2))+(2/3)*(2*j-3*loops)/loops; -R(j)/2 -R(j)/2 -R(j)*(1+1/sqrt(2))+(2/3)*(2*j-3*loops)/loops]; 
    hold off
    axes(ax)
    sn = surf(X*R(j),Y*R(j), R(j)*(Z+1)-(2/3)*(2*j-3*loops)/loops,0.*X);
    hold on
    ss = surf(R(j)*X,R(j)*Y,R(j)*(Z-1)+(2/3)*(2*j-3*loops)/loops,0.*X);
    scatter3(points(:,1),points(:,2),points(:,3),10*R(j),'x')
    plot3(points([1 2 4 3 1],1),points([1 2 4 3 1],2),points([1 2 4 3 1],3),'-.k')
    ylim([-1 1])
    xlim([-1 1])
    zlim([-2 2])
    sn.LineStyle = 'none';
    ss.LineStyle = 'none';
    sn.FaceAlpha = 0.3;
    ss.FaceAlpha = 0.3;
    view([30 40])
    daspect([1 1 1]);
    drawnow
    M(j) = getframe(gcf);
    writeVideo(v,M(j));
end
close(v);
%% 
h = figure;
h.Color = [.2 .2 .2];
points = 2E4;
r = 0.01;
start_rad = (r/10)*randn(points,1)+r;
start_angles = 2*pi*rand(points,2);
particle_0 = [start_rad(:).*sin(start_angles(:,1)).*cos(start_angles(:,2)) start_rad(:).*sin(start_angles(:,1)).*sin(start_angles(:,2)) start_rad(:).*cos(start_angles(:,1))];

north_angles = 2*pi*rand(points/2,2);
south_angles = 2*pi*rand(points/2,2);
north_momenta = [sin(north_angles(:,1)).*cos(north_angles(:,2)) sin(north_angles(:,1)).*sin(north_angles(:,2)) 1+cos(north_angles(:,1))];
south_momenta = [sin(south_angles(:,1)).*cos(south_angles(:,2)) sin(south_angles(:,1)).*sin(south_angles(:,2)) -1+cos(south_angles(:,1))];

point_loc = [0 r/sqrt(2) r/sqrt(2); 0 -r/sqrt(2) -r/sqrt(2); 0 r/sqrt(2) r/sqrt(2); 0 -r/sqrt(2) -r/sqrt(2)];
point_momenta = [0 1/sqrt(2) 1/sqrt(2)+1; 0 -1/sqrt(2) 1-1/sqrt(2); 0 1/sqrt(2) 1/sqrt(2)-1; 0 -1/sqrt(2) -1/sqrt(2)-1];
point_t = @(t) point_loc + t.*point_momenta;

ax1 = axes;
set(gca,'nextplot','replacechildren');
%h.Visible = 'off'
v = VideoWriter('GrowthSphere.avi');
v.FrameRate = 60;
open(v);


colour = [1.*ones(points/2,1).*(north_momenta(:,2)+1);-1.*ones(points/2,1).*(south_momenta(:,2)+1)];
colormap([winter(64)]);
momenta = [north_momenta ; south_momenta];
particle_t = @(t) particle_0 + t.*momenta;

tau = 0.01;
for t = 0:tau/239:tau
    part = particle_t(t);
    point = point_t(t);
    axes(ax1)
    hold off
    scatter3(ax1,part(:,1),part(:,2),part(:,3),0.5,sign(colour(:)))
    hold on
    plot3(point([1 2 4 3 1],1),point([1 2 4 3 1],2),point([1 2 4 3 1],3),'-.w')
    
    ax1.Color = [.2 .2 .2];
    ax1.XColor = [1 1 1];
    ax1.YColor = [1 1 1];
    ax1.ZColor = [1 1 1];
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.ZGrid = 'on';
    ax1.GridColor = [1 1 1];
    view([90,30])
    text(0,1.2*(point(1,2)),1.2*point(1,3),'$a^+$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(2,2)),1.2*point(2,3),'$a^-$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(3,2)),1.2*point(3,3),'$b^+$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(4,2)),1.2*point(4,3),'$b^-$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    xlim([-1.5*point(1,2) 1.5*point(1,2)])
    ylim([-1.5*point(1,2) 1.5*point(1,2)])
    zlim([-1.5*point(1,3) 1.5*point(1,3)])
    
    daspect([1 1 1])
    drawnow
    hold = getframe(gcf);
    writeVideo(v,hold);
end

particle_1 = particle_t(tau); point_1 = point_t(tau);
momenta(1:points/2,3) = momenta(1:points/2,3)-2; momenta(points/2+1:points,3) = momenta(points/2+1:points,3)+2;
point_momenta = [0 1/sqrt(2) 1/sqrt(2)-1; 0 -1/sqrt(2) -1-1/sqrt(2); 0 1/sqrt(2) 1/sqrt(2)+1; 0 -1/sqrt(2) -1/sqrt(2)+1];

point_t2 = @(t) point_1 + t.*point_momenta;
particle_t2 = @(t) particle_1 + t.*momenta;

tau = tau;
for t=0:tau/239:tau
    part = particle_t2(t);
    point = point_t2(t);
    hold off
    scatter3(ax1,part(:,1),part(:,2),part(:,3),0.5,sign(colour(:)))
    hold on
    plot3(point([1 2 4 3 1],1),point([1 2 4 3 1],2),point([1 2 4 3 1],3),'-.w')
    
    ax1.Color = [.2 .2 .2];
    ax1.XColor = [1 1 1];
    ax1.YColor = [1 1 1];
    ax1.ZColor = [1 1 1];
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.ZGrid = 'on';
    ax1.GridColor = [1 1 1];
    view([90,30])
    text(0,1.2*(point(1,2)),1.2*point(1,3),'$a^+$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(2,2)),1.2*point(2,3),'$a^-$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(3,2)),1.2*point(3,3),'$b^+$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(4,2)),1.2*point(4,3),'$b^-$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    xlim([-1.5*point(1,2) 1.5*point(1,2)])
    ylim([-1.5*point(1,2) 1.5*point(1,2)])
    zlim([-1.5*point(1,3) 1.5*point(1,3)])
    
    daspect([1 1 1])
    drawnow
    hold = getframe(gcf);
    writeVideo(v,hold);
end

particle_2 = particle_t2(tau);
point_2 = point_t2(tau);
colour = -colour;
particle_t3 = @(t) particle_2 + (t^3).*momenta;
point_t3 = @(t) point_2 + (t^3).*point_momenta;
colormap([autumn(32);cool(32)])
tau2 = 1.5;
for t=0:tau2/479:tau2
    part = particle_t3(t);
    point = point_t3(t);
    hold off
    scatter3(ax1,part(:,1),part(:,2),part(:,3),0.5,colour(:))
    hold on
    plot3(point([1 2 4 3 1],1),point([1 2 4 3 1],2),point([1 2 4 3 1],3),'-.w')
    ax1.Color = [.2 .2 .2];
    ax1.XColor = [1 1 1];
    ax1.YColor = [1 1 1];
    ax1.ZColor = [1 1 1];
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.ZGrid = 'on';
    ax1.GridColor = [1 1 1];
    xlim([-1.5*point(3,2) 1.5*point(3,2)])
    ylim([-1.5*point(3,2) 1.5*point(3,2)])
    zlim([-1.5*point(3,3) 1.5*point(3,3)])
    
    daspect([1 1 1])
    text(0,1.2*(point(1,2)),1.2*point(1,3),'$W$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(2,2)),1.2*point(2,3),'$Z$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(3,2)),1.2*point(3,3),'$Y$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    text(0,1.2*(point(4,2)),1.2*point(4,3),'$X$','Interpreter','latex','Color','w','FontSize',17,'HorizontalAlignment','center','VerticalAlignment','middle')
    
    view([90,30])

    drawnow
    hold = getframe(gcf);
    writeVideo(v,hold);
end
loops = 720;
for j = 1:loops
    view(90+j/2,15*(cos(2*pi*j/loops)+1))
    drawnow
    hold = getframe(gcf);
    writeVideo(v,hold);
end
close(v);
