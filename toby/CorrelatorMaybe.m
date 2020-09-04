tic
dYX = norm(Ycoord-Xcoord); dYW = norm(Ycoord-Wcoord); dYZ = norm(Ycoord-Zcoord);
dXW = norm(Xcoord-Wcoord); dXZ = norm(Xcoord-Zcoord); dWZ = norm(Wcoord-Zcoord);

dlist = -1*ones(1,cycles*Tmax);
index = 1;
for i = 1:cycles
    for j = 1:detectedYs(i)*detectedXs(i)
        dlist(index) = dYX;
        index = index+1;
    end
    for j=1:0.5*detectedYs(i)^2
        dlist(index) = 0;
        index = index+1;
    end
    for j=1:0.5*detectedXs(i)^2
        dlist(index) = 0;
        index = index+1;
    end
    for j=1:0.5*detectedWs(i)^2
        dlist(index) = 0;
        index = index+1;
    end
    for j=1:0.5*detectedZs(i)^2
        dlist(index) = 0;
        index = index+1;
    end
    for j=1:detectedYs(i)*detectedWs(i)
        dlist(index) = dYW;
        index = index+1;
    end
    for j = 1:detectedYs(i)*detectedZs(i)
        dlist(index) = dYZ;
        index = index+1;
    end
    for j=1:detectedXs(i)*detectedWs(i)
        dlist(index) = dXW;
        index = index+1;
    end
    for j = 1:detectedXs(i)*detectedZs(i)
        dlist(index) = dXZ;
        index = index+1;
    end
    for j=1:detectedWs(i)*detectedZs(i)
        dlist(index) = dWZ;
        index = index+1;
    end
end
dlist = dlist+1;
dlist = nonzeros(dlist');
dlist = dlist-1;


d2list = zeros(length(coordlist)^2,1);
index = 1;
for i = 1:length(coordlist)
    for j=1:length(coordlist)
        d2list(index) = norm(coordlist(i,:)-coordlist(j,:));
        index =index+1;
    end
end
d2list = nonzeros(d2list');
subplot(2,1,1)
histogram(dlist)
subplot(2,1,2)
histogram(d2list)
toc