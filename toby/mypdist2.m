function [YX,YZ,XW,WZ] = mypdist2(X)
    [n,~] = size(X);
    
    if n < 2
        YX = zeros(1,0);
        YZ = YX; XW = YX; WZ = YX;
        return;
    end
    
    % Divide the particles into their port location
    north = X(X(:,3)>0,:);
    south = X(X(:,3)<0,:);
    northeast = north(north(:,2)>0,:);
    northwest = north(north(:,2)<0,:);
    southeast = south(south(:,2)>0,:);
    southwest = south(south(:,2)<0,:);
    
    % Count the number of particles in each port region.
    [Nn,~] = size(north);
    [Sn,~] = size(south);
    [NEn,~] = size(northeast);
    [NWn,~] = size(northwest);
    [SEn,~] = size(southeast);
    [SWn,~] = size(southwest);
    
    % Create distance lists for each port-pair
    YX = zeros(1,NEn*NWn);
    YZ = zeros(1,NEn*Sn);
    XW = zeros(1,NWn*Sn);
    WZ = zeros(1,SEn*SWn);
    
    
    % Find the Y-Z pair back-to-back distances.
    k = 1;
    for i = 1:NEn
        dsq = (-northeast(i,1) - south(1:Sn,1)).^2 + (-northeast(i,2) - south(1:Sn,2)).^2 + (-northeast(i,3)-south(1:Sn,3)).^2;
        YZ(k:(k+Sn-1)) = sqrt(dsq);
        k = k + Sn;
    end
    
    % Find the X-W pair back-to-back distances.
    k =1;
    for i = 1:NWn
        dsq = (-northwest(i,1) - south(1:Sn,1)).^2 + (-northwest(i,2) - south(1:Sn,2)).^2 + (-northwest(i,3)-south(1:Sn,3)).^2;
        XW(k:(k+Sn-1)) = sqrt(dsq);
        k = k + Sn;
    end
    
    % Find the Y-X pair back-to-back distances.
    k = 1;
    for i = 1:NEn
        dsq = (-northeast(i,1) - northwest(1:NWn,1)).^2 + (-northeast(i,2) - northwest(1:NWn,2)).^2 + (2-northeast(i,3)-northwest(1:NWn,3)).^2;
        YX(k:(k+NWn-1)) = sqrt(dsq);
        k = k + NWn;
    end
    k=1;

    % Find the W-Z pair back-to-back distances.
    for i = 1:SEn
        dsq = (-southeast(i,1) - southwest(1:SWn,1)).^2 + (-southeast(i,2) - southwest(1:SWn,2)).^2 + (-2-southeast(i,3)-southwest(1:SWn,3)).^2;
        WZ(k:(k+SWn-1)) = sqrt(dsq);
        k = k + SWn;
    end
end