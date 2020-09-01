function colist_test = coordinatelist(lambda,phi,Tmax,cycles)
    
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

    colist_test = cell(1,cycles);

    % Blur factor:
    sigma = 0.05;
    for i = 1:cycles
        cordlist_temp = [];
        for j=1:detectedYs(i)
            cordlist_temp = [cordlist_temp; Ycoord+sigma*randn(1,3)];
        end
        for j=1:detectedXs(i)
            cordlist_temp = [cordlist_temp; Xcoord+ sigma*randn(1,3)];
        end
        for j=1:detectedWs(i)
            cordlist_temp = [cordlist_temp; Wcoord+ sigma*randn(1,3)];
        end
        for j=1:detectedZs(i)
                cordlist_temp = [cordlist_temp; Zcoord+ sigma*randn(1,3)];
        end
        colist_test{i} = cordlist_temp;
    end
end