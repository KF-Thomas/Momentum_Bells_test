% data = [0 0.92; pi/8 1.10; pi/2 1; 3*pi/4 0.94; pi 1.12; 3*pi/2 1.11];
data = [0 1.4 0.2; pi/4 1.5 0.2; pi/2 1.1 0.1; pi 1.7 0.1; 3*pi/2 0.45 0.1; 2*pi 1.45 0.1; 5*pi/4 0.8 0.1; 7*pi/4 1 0.1; 3*pi/4 1.5 0.1];
x = data(:,1);
y = data(:,2);
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz(:) .* circshift(yz(:),[1 0]) <= 0);     % Find zero-crossings
per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
xp = linspace(min(x),max(x));
fit = @(b,x)  b(1).*cos(x.*b(2) + 2*pi/b(3)).*(cos(x + 2*pi/b(3)).^2) + b(4);    % Function to fit
best_fit = fitnlm(x,y,fit,[0.8,0.02,-1.6,0.7625]);


[ysamp_val,ysamp_ci]=predict(best_fit,xp','Prediction','curve','Alpha',1-erf(1/sqrt(2))); %'Prediction','observation'



stfig('responce against phase');
clf
hold on
plot(xp,ysamp_val,'r','LineWidth',1.5)
drawnow
yl=ylim*1.1;
plot(xp,ysamp_ci,'color',[1,1,1].*0.5)
errorbar(x,y,data(:,3),'kx')
% plot(xp,fit(s,xp), 'r')
xlabel('$\phi$')
grid