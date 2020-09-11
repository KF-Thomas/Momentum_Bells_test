%% Generate our Bell Parameter off the correlations:
%% Lists:
phi = 0:pi/180:pi;
slist = 0.*phi;
%% Parameters:
cycles = 1E5; % C = number of test cycles simulated each time.
loss_rate = 8;% R = retention rate, or quantum efficiency.
dark_rate = 0;%  D = dark rate
blur = 0.05; % B = decay width of gaussian blur.
lambda = 0.1; % L = mode occupancy.
%% Evaluation for each phase:
for p = 1:length(phi)
    colist_test = coordinatelist(lambda,phi(p),8,cycles,loss_rate,dark_rate,blur);
    dlist = [];
    for i = 1:cycles
        temp = colist_test{i};
        dlist = [dlist pdist(temp)];
    end
    colist2 = [];
    for i = 1:cycles
        colist2 = [colist2; colist_test{i}];
    end
    d2list = pdist(colist2);
    
    edges = 0:.01:3.2;
    
    h1 = histogram(dlist,edges,'Normalization','probability');
    singlerun = h1.Values;
    
    h2 = histogram(d2list,edges,'Normalization','probability');
    totalrun = h2.Values;
    hlist = singlerun ./ totalrun;
    
    % Average over nearby values to smoothen out curve or take the values
    % only at exactly the edge / exactly the diagonal:
    edge = mean(hlist(196:206));
    diag = mean(hlist(278:288));
    %edge = hlist(201);
    %diag = hlist(283);
    
    S = (2*edge-2*diag)/(2*edge+2*diag);
    slist(p) = S;
end
%% Generate Plot / Image:
plot(phi,slist)
xticks([0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3,2*pi])
xticklabels(["0","\pi/3","2\pi/3","\pi","4\pi/3","5\pi/3","2\pi"])
xlabel("Phase ($\phi$)","Interpreter",'latex')
ylabel("Bell Parameter ($E(\phi)$)","Interpreter","latex")
title("Bell Parameter vs Phase,$\lambda = 0.1$, $R = 8$\%, $D = 0$, B = 0.05, $C = 10^5$","Interpreter","latex")

%% Write to File:
A = [phi; slist];
data = fopen('R8D00B005C1e5L01.txt','w');
fprintf(data,'%6s %12s\n','\phi','E(\phi)');
fprintf(data,'%6.4f %12.8f\n',A);
fclose(data);