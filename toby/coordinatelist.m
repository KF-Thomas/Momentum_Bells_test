% Returns a cell array with each cell containing the list of detected
% particle coordinates for each simulated experimental run.
%% Function:
function colist_test = coordinatelist(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur)
    % Call the divided probability list from problist.m
    [p,event,~] = problist(lambda,phi,Tmax);
    % Lists of each event for each trial:
    ys = zeros(cycles,1); xs = ys; ws = ys; zs = ys;
    % Find what event each trial corresponds to and note where these events
    % correlate to detections:
    for runs = 1:cycles
        n = rand;
        i = 1;
        while true
            if n<p(i+1)
                ys(runs) = event(i,1);
                xs(runs) = event(i,2);
                ws(runs) = event(i,3);
                zs(runs) = event(i,4);
                break
            end
            i=i+1;
        end
    end
    %% Quantum Efficiency:
    % Cull data, keeping only the efficiency rating:
    thresh = quantum_efficiency +.01;
    keeplist = randi(100,cycles,1)<thresh;
    survivedYs = ys.*keeplist;
    survivedXs = xs.*keeplist;
    survivedWs = ws.*keeplist;
    survivedZs = zs.*keeplist;

    %% Dark Counts:
    % Scale percentage by 10^6 so that small percentages can affect it.
    thresh = dark_rate*1E8+.1;
    darkYs = randi(1E8,cycles,1)<thresh;
    darkXs = randi(1E8,cycles,1)<thresh;
    darkWs = randi(1E8,cycles,1)<thresh;
    darkZs = randi(1E8,cycles,1)<thresh;
    detectedYs = survivedYs + darkYs;
    detectedXs = survivedXs + darkXs;
    detectedWs = survivedWs + darkWs;
    detectedZs = survivedZs + darkZs;

    %% Place detected particles at a position:
    Ycoord = [1 1 0];
    Xcoord = [-1 1 0];
    Wcoord = [1 -1 0];
    Zcoord = [-1 -1 0];

    colist_test = cell(1,cycles);

    % Blur factor:
    sigma = blur;
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