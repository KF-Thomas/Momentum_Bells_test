phi = 0:pi/180:pi;
slist = 0.*phi;
cycles = 2E4;
for p = 1:length(phi)
    colist_test = coordinatelist(0.8,phi(p),8,cycles);
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
    %figure(1);
    %subplot(2,1,1)
    h1 = histogram(dlist,edges,'Normalization','probability');
    singlerun = h1.Values;
    %subplot(2,1,2)
    h2 = histogram(d2list,edges,'Normalization','probability');
    totalrun = h2.Values;
    %figure(2);
    hlist = singlerun ./ totalrun;
    %plot(edges(1:end-1),hlist)
    edge = mean(hlist(195:205));
    diag = mean(hlist(278:288));
    S = (2*edge-2*diag)/(2*edge+2*diag);
    slist(p) = S;
end
plot(phi,slist)
xticks([0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3,2*pi])
xticklabels({"0","\pi/3","2\pi/3","\pi","4\pi/3","5\pi/3","2\pi"})
