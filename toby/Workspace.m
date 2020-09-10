edges = 0:.01:3.3;
h1 = histogram(dlist,edges,'Normalization','probability');
singlerun = h1.Values;
h2 = histogram(d2list,edges,'Normalization','probability');
Totalrun = h2.Values;

hlist = singlerun ./ Totalrun;
plot(hlist)