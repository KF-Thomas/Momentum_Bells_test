%Generate our Bell Parameter off the correlations:
%Lists:
phi = 0:pi/180:pi;
slist = 0.*phi;

%Parameters:
cycles = 1E6; % C = number of test cycles simulated each time.
loss_rate = 100; % R = retention rate, or quantum efficiency.
dark_rate = 1E-3; % D = dark rate
blur = 0.05; % B = decay width of gaussian blur.
lambda = 0.1; % L = mode occupancy.

phi = pi/2;

R = 1;
theta = pi/4;

Ycoord = [R*cos(theta) R*(1+sin(theta)) 0];
Xcoord = [R*cos(theta) R*(sin(theta)-1) 0];
Wcoord = [-R*cos(theta) R*(1-sin(theta)) 0];
Zcoord = [-R*cos(theta) -R*(1+sin(theta)) 0];

%Evaluation for each phase:
%for p = 1:length(phi)
colist_test = coordinatelist(lambda,phi,10,cycles,loss_rate,dark_rate,blur,Ycoord,Xcoord,Wcoord,Zcoord);
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

edges = 0:.01:4;
    
histogram(dlist,edges,'Normalization','probability');
return
singlerun = h1.Values;
return
h2 = histogram(d2list,edges,'Normalization','probability');
totalrun = h2.Values;
hlist = singlerun ./ totalrun;

% Average over nearby values to smoothen out curve or take the values
% only at exactly the edge / exactly the diagonal:
%side = 
edge = mean(hlist(196:206));
diag = mean(hlist(278:288));
%edge = hlist(201);
%diag = hlist(283);

S = (2*edge-2*diag)/(2*edge+2*diag);
slist(p) = S;
%end
