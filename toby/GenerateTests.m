phi = pi/2;
Tmax = 15;
lambda = 0.8;
cycles = 1E4;

[p,event,cap] = problist(lambda,phi,Tmax);

m = zeros(cycles,1); ys = m; xs = m; ws = m; zs = m;

for runs = 1:cycles
    n = rand;
    i = 1;
    while true
        if n<p(i+1)
            m(runs) = i;
            ys(runs) = event(i,1);
            xs(runs) = event(i,2);
            ws(runs) = event(i,3);
            zs(runs) = event(i,4);
            break
        end
        i=i+1;
    end
end
% Lost Counts:
lostlist = randi(100,cycles,1)<8.01;
survivedm = m.*lostlist;
survivedYs = ys.*lostlist;
survivedXs = xs.*lostlist;
survivedWs = ws.*lostlist;
survivedZs = zs.*lostlist;

% Dark Counts:
darkYs = randi(1000,cycles,1)<1.01;
darkXs = randi(1000,cycles,1)<1.01;
darkWs = randi(1000,cycles,1)<1.01;
darkZs = randi(1000,cycles,1)<1.01;
detectedYs = survivedYs + darkYs;
detectedXs = survivedXs + darkXs;
detectedWs = survivedWs + darkWs;
detectedZs = survivedZs + darkZs;



coordlist = zeros(sum(detectedYs)+sum(detectedXs)+sum(detectedWs)+sum(detectedZs),3);
pos = 1;
Ycoord = [1 1 0];
Xcoord = [-1 1 0];
Wcoord = [1 -1 0];
Zcoord = [-1 -1 0];

% Blur factor:
sigma = 0.05;

for j = 1:sum(detectedYs)
    coordlist(pos,:) = Ycoord + sigma*randn(1,3);
    pos = pos+1;
end
for j = 1:sum(detectedXs)
    coordlist(pos,:) = Wcoord + sigma*randn(1,3);
    pos = pos+1;
end 

for j = 1:sum(detectedWs)
    coordlist(pos,:) = Xcoord+sigma*randn(1,3);
    pos = pos+1;
end 

for j = 1:sum(detectedZs)
    coordlist(pos,:) = Zcoord+sigma*randn(1,3);
    pos = pos+1;
end 

% Show the detections with the blur factor.
scatter3(coordlist(:,1),coordlist(:,2),coordlist(:,3))
zlim([-1 1])