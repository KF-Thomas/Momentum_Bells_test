function [YX,YZ,XW,WZ] = mypdist2(X)
    [n,~] = size(X);
    
    if n < 2
        YX = zeros(1,0);
        YZ = YX; XW = YX; WZ = YX;
        return;
    end
    
    north = X(X(:,3)>0,:);
    south = X(X(:,3)<0,:);
    northeast = north(north(:,2)>0,:);
    northwest = north(north(:,2)<0,:);
    southeast = south(south(:,2)>0,:);
    southwest = south(south(:,2)<0,:);
    
    [Nn,~] = size(north);
    [Sn,~] = size(south);
    [NEn,~] = size(northeast);
    [NWn,~] = size(northwest);
    [SEn,~] = size(southeast);
    [SWn,~] = size(southwest);
    
    YX = zeros(1,Nn*(Nn-1)/2);
    YZ = zeros(1,NEn*Sn);
    XW = zeros(1,NWn*Sn);
    WZ = zeros(1,Sn*(Sn-1)/2);
    
    k = 1;
    for i = 1:NEn
        dsq = (-northeast(i,1) - south(1:Sn,1)).^2 + (-northeast(i,2) - south(1:Sn,2)).^2 + (-northeast(i,3)-south(1:Sn,3)).^2;
        YZ(k:(k+Sn-1)) = sqrt(dsq);
        k = k + Sn;
    end
    k =1;
    for i = 1:NWn
        dsq = (-northwest(i,1) - south(1:Sn,1)).^2 + (-northwest(i,2) - south(1:Sn,2)).^2 + (-northwest(i,3)-south(1:Sn,3)).^2;
        XW(k:(k+Sn-1)) = sqrt(dsq);
        k = k + Sn;
    end
    k = 1;
    for i = 1:NEn
        dsq = (-northeast(i,1) - northwest(1:NWn,1)).^2 + (-northeast(i,2) - northwest(1:NWn,2)).^2 + (2-northeast(i,3)-northwest(1:NWn,3)).^2;
        YX(k:(k+NWn-1)) = sqrt(dsq);
        k = k + NWn;
    end
    k=1;
%    for i = 1:Nn-1
%        dsq = (-north(i,1)-north((i+1):Nn,1)).^2 +(-north(i,2)-north((i+1):Nn,2)).^2 + (2-north(i,3)-north((i+1):Nn,3)).^2;
%        YX(k:(k+Nn-i-1)) = sqrt(dsq);
%        k = k+Nn-i;
%    end
   for i = 1:SEn
        dsq = (-southeast(i,1) - southwest(1:SWn,1)).^2 + (-southeast(i,2) - southwest(1:SWn,2)).^2 + (-2-southeast(i,3)-southwest(1:SWn,3)).^2;
        WZ(k:(k+SWn-1)) = sqrt(dsq);
        k = k + SWn;
    end
 
    k=1;
    %for i = 1:Sn-1
    %    dsq = (-south(i,1)-south((i+1):Sn,1)).^2 +(-south(i,2)-south((i+1):Sn,2)).^2 + (-2-south(i,3)-south((i+1):Sn,3)).^2;
    %    WZ(k:(k+Sn-i-1)) = sqrt(dsq);
    %    k=k+Sn-i;
    %end

end